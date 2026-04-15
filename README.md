# BioAuto

BioAuto is a small local web app for lightweight bioinformatics automation and open-source inhibitor lookup.

## Gene inhibitor search

The main screen lets you type a human gene symbol or gene-plus-variant, such as `PARP1`, `EGFR L858R`, `BRAF V600E`, or `ESR1`, and returns up to 50 candidate inhibitors with mechanism evidence and citation links.

Strict matching rules:

- The first token must exactly match a human ChEMBL target gene symbol. Fuzzy spellings such as `parpy1` will not be silently mapped to `PARP1`.
- If you enter a disease or tissue context, such as `fallopian tube cancer`, exact-context mode is enabled by default.
- In exact-context mode, an inhibitor is shown only when an open-access Europe PMC citation explicitly matches the typed context together with the gene and inhibitor.
- The app will return no contextual hit rather than substituting related diseases such as ovarian cancer for fallopian tube cancer.
- Use the exclusion field for hard separation, for example context `fallopian tube cancer` plus exclude `ovarian`.

Sources used by the app:

- ChEMBL open APIs for target matching, curated drug mechanisms, molecule metadata, and mechanism references.
- PubMed E-utilities for PubMed citation titles attached to ChEMBL mechanism records.
- Europe PMC open-access search for supplemental open-paper links when additional citation context is available.

Important limitations:

- This is a research-assist tool, not medical advice or a replacement for manual review.
- A result means an open database or open literature source connects the molecule to an inhibitory/antagonist/degrader mechanism for the exact matched target.
- Always verify the primary paper, target isoform, species, assay system, and disease context before using a result in an experiment or report.
- Some genes, including tumor suppressors such as `BRCA1`, may not have direct small-molecule inhibitors in ChEMBL. In those cases, search pathway or receptor targets such as `ESR1`, `PARP1`, or `EGFR` depending on the biology.

## Included sequence workflows

- FASTA summary: built into the app, no external tool required.
- seqkit stats: requires `seqkit` on `PATH`.
- FastQC quality check: requires `fastqc` on `PATH`.
- BLASTN search: requires `blastn` on `PATH` and a local BLAST database path.

## Run the app

```powershell
python app.py
```

Then open:

```text
http://127.0.0.1:8000
```

Try the built-in sequence workflow with:

```text
sample-data/example.fasta
```

## Add another command-line tool

Edit `workflows.json` and add a workflow entry. Use these placeholders in commands:

- `{input}`: uploaded input file
- `{output_dir}`: folder for this run's outputs
- `{threads}`: thread count from the form
- `{database}`: BLAST database path from the form

Example:

```json
{
  "id": "my_tool",
  "name": "My tool",
  "description": "Runs my command-line analysis.",
  "executable": "mytool",
  "command": ["mytool", "--input", "{input}", "--out", "{output_dir}/results.txt"],
  "input_hint": "FASTA"
}
```

Run outputs are written under `runs/`.
