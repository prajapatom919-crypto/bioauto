from __future__ import annotations

import cgi
import html
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parent
STATIC_DIR = ROOT / "static"
RUNS_DIR = ROOT / "runs"
WORKFLOWS_FILE = ROOT / "workflows.json"
CH_EMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
EUROPE_PMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
REQUEST_TIMEOUT = 18
MAX_INHIBITORS = 50
HTTP_CACHE: dict[str, dict] = {}




def fetch_json(url: str, params: dict | None = None, timeout: int = REQUEST_TIMEOUT) -> dict:
    if params:
        query = urllib.parse.urlencode({key: value for key, value in params.items() if value is not None})
        url = f"{url}?{query}"

    if url in HTTP_CACHE:
        return HTTP_CACHE[url]

    request = urllib.request.Request(url, headers={"User-Agent": "BioAuto inhibitor search/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Open API request failed with HTTP {exc.code}: {url}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not reach an open bioinformatics API: {exc.reason}") from exc

    HTTP_CACHE[url] = payload
    return payload


def parse_gene_query(value: str) -> tuple[str, str]:
    cleaned = re.sub(r"[^A-Za-z0-9_.:>\-\s]+", " ", value.strip())
    parts = cleaned.split()
    if not parts:
        raise ValueError("Type a gene symbol or gene variant, for example ESR1, EGFR L858R, or BRAF V600E.")
    gene = re.sub(r"[^A-Za-z0-9_.-]+", "", parts[0]).upper()
    variant = " ".join(parts[1:]).upper()
    if not gene:
        raise ValueError("Type an official gene symbol first, then an optional variant, for example EGFR L858R.")
    if len(gene) > 40:
        raise ValueError("Gene symbol is too long. Use an official short gene symbol such as ESR1, EGFR, or PARP1.")
    return gene, variant


def normalize_context(value: str | None) -> str:
    context = re.sub(r"\s+", " ", (value or "").strip())
    if len(context) > 120:
        raise ValueError("Disease or tissue context is too long. Use a concise phrase such as fallopian tube cancer.")
    return context


def normalize_exclusions(value: str | None) -> list[str]:
    raw_terms = re.split(r"[,;]", value or "")
    terms = []
    for term in raw_terms:
        cleaned = re.sub(r"\s+", " ", term.strip())
        if cleaned:
            if len(cleaned) > 60:
                raise ValueError("Each exclusion term must be concise, for example ovarian.")
            terms.append(cleaned)
    return terms[:8]


def component_symbols(target: dict) -> set[str]:
    symbols: set[str] = set()
    for component in target.get("target_components") or []:
        for synonym in component.get("target_component_synonyms") or []:
            if synonym.get("syn_type", "").startswith("GENE_SYMBOL"):
                symbols.add(str(synonym.get("component_synonym", "")).upper())
    return symbols


def find_chembl_target(gene: str) -> dict:
    payload = fetch_json(f"{CH_EMBL_BASE}/target/search.json", {"q": gene, "limit": 50})
    targets = payload.get("targets") or []
    exact_human = [
        target
        for target in targets
        if target.get("organism") == "Homo sapiens" and gene in component_symbols(target)
    ]

    def score(target: dict) -> tuple[int, int]:
        single_protein = target.get("target_type") == "SINGLE PROTEIN"
        return (int(single_protein), int(target.get("species_group_flag") is False))

    if not exact_human:
        raise ValueError(
            f"No exact human ChEMBL target matched gene symbol {gene}. "
            "Use the official gene symbol first, for example PARP1 rather than parpy1."
        )
    return sorted(exact_human, key=score, reverse=True)[0]


def is_inhibitory_mechanism(mechanism: dict) -> bool:
    action = str(mechanism.get("action_type") or "").upper()
    text = " ".join(
        str(mechanism.get(field) or "")
        for field in ("mechanism_of_action", "mechanism_comment", "selectivity_comment")
    ).upper()
    inhibitory_actions = {
        "ANTAGONIST",
        "INHIBITOR",
        "BLOCKER",
        "INVERSE AGONIST",
        "NEGATIVE ALLOSTERIC MODULATOR",
        "DEGRADER",
    }
    inhibitory_words = ("ANTAGONIST", "INHIBITOR", "INHIBITS", "BLOCKER", "DEGRADER", "SUPPRESS")
    return action in inhibitory_actions or any(word in text for word in inhibitory_words)


def fetch_mechanisms(target_chembl_id: str) -> list[dict]:
    payload = fetch_json(
        f"{CH_EMBL_BASE}/mechanism.json",
        {"target_chembl_id": target_chembl_id, "limit": 1000},
    )
    return [item for item in payload.get("mechanisms") or [] if is_inhibitory_mechanism(item)]


def fetch_molecule(molecule_chembl_id: str) -> dict:
    return fetch_json(f"{CH_EMBL_BASE}/molecule/{urllib.parse.quote(molecule_chembl_id)}.json")


def citation_title_from_pubmed(pmid: str) -> str | None:
    if not pmid.isdigit():
        return None
    payload = fetch_json(
        NCBI_ESUMMARY,
        {"db": "pubmed", "id": pmid, "retmode": "json"},
        timeout=12,
    )
    return (((payload.get("result") or {}).get(pmid) or {}).get("title") or None)


def open_article_hits(
    gene: str,
    inhibitor_name: str,
    limit: int = 3,
    variant: str = "",
    context: str = "",
    exclude_terms: list[str] | None = None,
) -> list[dict]:
    query_parts = ["OPEN_ACCESS:y", f"({gene})", f'("{inhibitor_name}" OR inhibitor OR antagonist)']
    if variant:
        query_parts.append(f'"{variant}"')
    if context:
        query_parts.append(f'"{context}"')
    for term in exclude_terms or []:
        query_parts.append(f'NOT "{term}"')
    query = " AND ".join(query_parts)
    payload = fetch_json(
        EUROPE_PMC_BASE,
        {"query": query, "format": "json", "pageSize": limit, "resultType": "core"},
        timeout=14,
    )
    hits = []
    for item in ((payload.get("resultList") or {}).get("result") or [])[:limit]:
        pmid = item.get("pmid") or item.get("id")
        doi = item.get("doi")
        hits.append(
            {
                "title": item.get("title") or "Open-access article",
                "source": item.get("journalTitle") or item.get("source") or "Europe PMC",
                "year": item.get("pubYear"),
                "pmid": pmid,
                "doi": doi,
                "url": item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url")
                if item.get("fullTextUrlList")
                else (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else (f"https://doi.org/{doi}" if doi else None)),
                "evidence_source": "Europe PMC exact-context open-access search" if context else "Europe PMC open-access search",
            }
        )
    return hits


def format_mechanism_refs(refs: list[dict]) -> list[dict]:
    citations = []
    seen = set()
    for ref in refs or []:
        ref_type = ref.get("ref_type") or "Reference"
        ref_id = str(ref.get("ref_id") or "").strip()
        ref_url = ref.get("ref_url")
        key = (ref_type, ref_id, ref_url)
        if key in seen:
            continue
        seen.add(key)
        title = None
        if ref_type == "PubMed":
            try:
                title = citation_title_from_pubmed(ref_id)
            except Exception:
                title = None
        citations.append(
            {
                "title": title or f"{ref_type}: {ref_id}".strip(),
                "source": ref_type,
                "pmid": ref_id if ref_type == "PubMed" and ref_id.isdigit() else None,
                "url": ref_url or (f"https://pubmed.ncbi.nlm.nih.gov/{ref_id}/" if ref_type == "PubMed" and ref_id.isdigit() else None),
                "evidence_source": "ChEMBL mechanism reference",
            }
        )
    return citations


def inhibitor_search(
    gene_input: str,
    limit_input: str | None = None,
    context_input: str | None = None,
    require_context_input: str | None = None,
    exclude_input: str | None = None,
) -> dict:
    gene, variant = parse_gene_query(gene_input)
    context = normalize_context(context_input)
    exclude_terms = normalize_exclusions(exclude_input)
    require_context = context and str(require_context_input or "true").lower() != "false"
    try:
        limit = min(max(int(limit_input or MAX_INHIBITORS), 1), MAX_INHIBITORS)
    except ValueError:
        limit = MAX_INHIBITORS

    target = find_chembl_target(gene)
    mechanisms = fetch_mechanisms(target["target_chembl_id"])
    grouped: dict[str, dict] = {}

    for mechanism in mechanisms:
        molecule_id = mechanism.get("parent_molecule_chembl_id") or mechanism.get("molecule_chembl_id")
        if not molecule_id:
            continue
        entry = grouped.setdefault(
            molecule_id,
            {
                "molecule_chembl_id": molecule_id,
                "mechanisms": [],
                "citations": [],
                "max_phase": mechanism.get("max_phase"),
            },
        )
        entry["mechanisms"].append(
            {
                "action_type": mechanism.get("action_type"),
                "mechanism_of_action": mechanism.get("mechanism_of_action"),
                "comment": mechanism.get("mechanism_comment"),
                "direct_interaction": bool(mechanism.get("direct_interaction")),
            }
        )
        entry["citations"].extend(format_mechanism_refs(mechanism.get("mechanism_refs") or []))
        if mechanism.get("max_phase") is not None:
            entry["max_phase"] = max(entry.get("max_phase") or 0, mechanism.get("max_phase") or 0)

    results = []
    for molecule_id, entry in sorted(grouped.items(), key=lambda item: item[1].get("max_phase") or 0, reverse=True)[:limit]:
        molecule = {}
        try:
            molecule = fetch_molecule(molecule_id)
        except Exception:
            molecule = {"molecule_chembl_id": molecule_id}

        citations = []
        citation_keys = set()
        for citation in entry["citations"]:
            key = (citation.get("title"), citation.get("url"), citation.get("pmid"))
            if key not in citation_keys:
                citation_keys.add(key)
                citations.append(citation)

        name = molecule.get("pref_name") or molecule.get("molecule_pref_name") or molecule_id
        context_citations = []
        if context and name != molecule_id:
            try:
                context_citations = open_article_hits(gene, name, 5, variant=variant, context=context, exclude_terms=exclude_terms)
            except Exception:
                context_citations = []

        if require_context and not context_citations:
            continue

        if len(citations) < 3 and name != molecule_id and not require_context:
            try:
                for hit in open_article_hits(gene, name, 3 - len(citations), variant=variant, exclude_terms=exclude_terms):
                    key = (hit.get("title"), hit.get("url"), hit.get("pmid"))
                    if key not in citation_keys:
                        citation_keys.add(key)
                        citations.append(hit)
            except Exception:
                pass

        for hit in context_citations:
            key = (hit.get("title"), hit.get("url"), hit.get("pmid"))
            if key not in citation_keys:
                citation_keys.add(key)
                citations.append(hit)

        results.append(
            {
                "name": name,
                "chembl_id": molecule_id,
                "max_phase": molecule.get("max_phase") or entry.get("max_phase"),
                "first_approval": molecule.get("first_approval"),
                "molecule_type": molecule.get("molecule_type"),
                "mechanisms": entry["mechanisms"],
                "citations": citations,
                "context_citation_count": len(context_citations),
                "context_match_required": bool(require_context),
                "chembl_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{molecule_id}/",
            }
        )

    return {
        "ok": True,
        "gene": gene,
        "variant": variant,
        "context": context,
        "exclude_terms": exclude_terms,
        "require_context": bool(require_context),
        "target": {
            "chembl_id": target.get("target_chembl_id"),
            "name": target.get("pref_name"),
            "organism": target.get("organism"),
            "type": target.get("target_type"),
            "symbols": sorted(component_symbols(target)),
            "url": f"https://www.ebi.ac.uk/chembl/target_report_card/{target.get('target_chembl_id')}/",
        },
        "count": len(results),
        "max_results": limit,
        "sources": [
            "ChEMBL open target, molecule, and drug-mechanism APIs",
            "PubMed E-utilities for mechanism-reference titles",
            "Europe PMC open-access search for supplemental open literature links",
        ],
        "disclaimer": "Strict mode: exact human gene-symbol target matching is required. If a disease/tissue context is supplied, results are shown only when open-access literature explicitly matches that context unless exact-context mode is turned off.",
        "results": results,
    }

def load_workflows() -> dict:
    with WORKFLOWS_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def json_response(handler: SimpleHTTPRequestHandler, payload: dict, status: int = 200) -> None:
    body = json.dumps(payload, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def safe_filename(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    return cleaned or "input.txt"


def root_relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_name: str | None = None
    chunks: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, "".join(chunks)))
                current_name = line[1:].strip() or f"sequence_{len(records) + 1}"
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line).upper())

    if current_name is not None:
        records.append((current_name, "".join(chunks)))
    return records


def fasta_summary(input_path: Path, output_dir: Path) -> dict:
    records = parse_fasta(input_path)
    if not records:
        raise ValueError("No FASTA records were found. Upload a file with headers that start with '>'.")

    rows = []
    total_length = 0
    gc_count = 0
    n_count = 0

    for name, sequence in records:
        length = len(sequence)
        gc = sequence.count("G") + sequence.count("C")
        n_bases = sequence.count("N")
        total_length += length
        gc_count += gc
        n_count += n_bases
        rows.append(
            {
                "name": name,
                "length": length,
                "gc_percent": round((gc / length) * 100, 2) if length else 0,
                "n_bases": n_bases,
            }
        )

    summary = {
        "records": len(records),
        "total_length": total_length,
        "average_length": round(total_length / len(records), 2),
        "gc_percent": round((gc_count / total_length) * 100, 2) if total_length else 0,
        "n_bases": n_count,
        "sequences": rows,
    }

    report_path = output_dir / "fasta-summary.json"
    html_path = output_dir / "fasta-summary.html"
    report_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    html_path.write_text(render_fasta_report(summary), encoding="utf-8")

    return {
        "summary": summary,
        "stdout": json.dumps(summary, indent=2),
        "stderr": "",
        "files": [root_relative(report_path), root_relative(html_path)],
    }


def render_fasta_report(summary: dict) -> str:
    rows = "\n".join(
        "<tr>"
        f"<td>{html.escape(row['name'])}</td>"
        f"<td>{row['length']}</td>"
        f"<td>{row['gc_percent']}</td>"
        f"<td>{row['n_bases']}</td>"
        "</tr>"
        for row in summary["sequences"]
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>FASTA Summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #172026; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
    th, td {{ border: 1px solid #c9d4d8; padding: 8px; text-align: left; }}
    th {{ background: #eef5f1; }}
  </style>
</head>
<body>
  <h1>FASTA Summary</h1>
  <p><strong>Records:</strong> {summary['records']}</p>
  <p><strong>Total length:</strong> {summary['total_length']}</p>
  <p><strong>Average length:</strong> {summary['average_length']}</p>
  <p><strong>GC percent:</strong> {summary['gc_percent']}</p>
  <table>
    <thead><tr><th>Name</th><th>Length</th><th>GC %</th><th>N bases</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</body>
</html>
"""


def run_external_tool(workflow: dict, input_path: Path, output_dir: Path, parameters: dict) -> dict:
    executable = workflow["executable"]
    if shutil.which(executable) is None:
        raise FileNotFoundError(
            f"'{executable}' is not installed or is not on PATH. Install it, then run this workflow again."
        )

    context = {
        "input": str(input_path),
        "output_dir": str(output_dir),
        "threads": str(parameters.get("threads") or workflow.get("default_threads", 2)),
        "database": str(parameters.get("database") or ""),
    }

    command = [part.format(**context) for part in workflow["command"]]
    started = time.time()
    result = subprocess.run(command, capture_output=True, text=True, cwd=ROOT)
    duration = round(time.time() - started, 2)

    log_path = output_dir / "command-log.txt"
    log_path.write_text(
        "Command:\n"
        + " ".join(command)
        + "\n\nExit code:\n"
        + str(result.returncode)
        + "\n\nSTDOUT:\n"
        + result.stdout
        + "\n\nSTDERR:\n"
        + result.stderr,
        encoding="utf-8",
    )

    if result.returncode != 0:
        raise RuntimeError(f"{workflow['name']} failed. See command-log.txt in the run folder.")

    files = [root_relative(path) for path in output_dir.rglob("*") if path.is_file()]
    return {
        "summary": {"duration_seconds": duration, "command": command},
        "stdout": result.stdout,
        "stderr": result.stderr,
        "files": files,
    }


class BioAutoHandler(SimpleHTTPRequestHandler):
    def translate_path(self, path: str) -> str:
        parsed = urlparse(path)
        if parsed.path.startswith("/runs/"):
            return str(ROOT / parsed.path.lstrip("/"))
        if parsed.path in ("/", "/index.html"):
            return str(STATIC_DIR / "index.html")
        return str(STATIC_DIR / parsed.path.lstrip("/"))

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/api/inhibitors":
            query = urllib.parse.parse_qs(parsed.query)
            try:
                payload = inhibitor_search(
                    query.get("gene", [""])[0],
                    query.get("limit", [str(MAX_INHIBITORS)])[0],
                    query.get("context", [""])[0],
                    query.get("require_context", ["true"])[0],
                    query.get("exclude", [""])[0],
                )
                json_response(self, payload)
            except Exception as exc:
                json_response(self, {"ok": False, "error": str(exc)}, status=400)
            return

        if parsed.path == "/api/workflows":
            workflows = load_workflows()
            for workflow in workflows["workflows"]:
                executable = workflow.get("executable")
                workflow["available"] = executable == "builtin" or shutil.which(executable) is not None
            json_response(self, workflows)
            return
        super().do_GET()

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path != "/api/run":
            json_response(self, {"error": "Unknown endpoint"}, status=404)
            return

        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={
                "REQUEST_METHOD": "POST",
                "CONTENT_TYPE": self.headers.get("Content-Type", ""),
            },
        )

        workflow_id = form.getfirst("workflow", "")
        parameters = json.loads(form.getfirst("parameters", "{}") or "{}")
        upload = form["input_file"] if "input_file" in form else None
        if upload is None or not upload.filename:
            json_response(self, {"error": "Upload an input FASTA/FASTQ file."}, status=400)
            return

        workflows = {item["id"]: item for item in load_workflows()["workflows"]}
        if workflow_id not in workflows:
            json_response(self, {"error": "Choose a valid workflow."}, status=400)
            return

        run_id = time.strftime("%Y%m%d-%H%M%S")
        run_dir = RUNS_DIR / run_id
        output_dir = run_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = run_dir / safe_filename(upload.filename)

        with input_path.open("wb") as handle:
            handle.write(upload.file.read())

        workflow = workflows[workflow_id]
        try:
            if workflow["executable"] == "builtin":
                result = fasta_summary(input_path, output_dir)
            else:
                result = run_external_tool(workflow, input_path, output_dir, parameters)
        except Exception as exc:
            error_path = output_dir / "error.txt"
            error_path.write_text(str(exc), encoding="utf-8")
            json_response(
                self,
                {
                    "ok": False,
                    "run_id": run_id,
                    "error": str(exc),
                    "files": [root_relative(input_path), root_relative(error_path)],
                },
                status=500,
            )
            return

        json_response(
            self,
            {
                "ok": True,
                "run_id": run_id,
                "workflow": workflow["name"],
                "input": root_relative(input_path),
                **result,
            },
        )


def main() -> None:
    port = int(os.environ.get("PORT", "8000"))
    RUNS_DIR.mkdir(exist_ok=True)
    server = ThreadingHTTPServer(("0.0.0.0", port), BioAutoHandler)
    print(f"BioAuto is running at http://127.0.0.1:{port}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
        sys.exit(0)


if __name__ == "__main__":
    main()
