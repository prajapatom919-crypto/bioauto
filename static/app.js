const workflowSelect = document.querySelector("#workflow");
const details = document.querySelector("#workflow-details");
const databaseRow = document.querySelector("#database-row");
const form = document.querySelector("#run-form");
const button = document.querySelector("#run-button");
const statusBox = document.querySelector("#status");
const stdoutBox = document.querySelector("#stdout");
const filesBox = document.querySelector("#files");
const tabs = document.querySelectorAll(".tab");
const panels = document.querySelectorAll(".panel");
const inhibitorForm = document.querySelector("#inhibitor-form");
const inhibitorButton = document.querySelector("#inhibitor-button");
const inhibitorStatus = document.querySelector("#inhibitor-status");
const inhibitorCount = document.querySelector("#inhibitor-count");
const targetCard = document.querySelector("#target-card");
const inhibitorList = document.querySelector("#inhibitor-list");

let workflows = [];

function selectedWorkflow() {
  return workflows.find((item) => item.id === workflowSelect.value);
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function renderWorkflow() {
  const workflow = selectedWorkflow();
  if (!workflow) return;

  databaseRow.classList.toggle("hidden", !workflow.requires_database);
  details.innerHTML = `
    <strong>${escapeHtml(workflow.name)}</strong><br>
    ${escapeHtml(workflow.description)}<br>
    <span>Input: ${escapeHtml(workflow.input_hint)}</span><br>
    ${
      workflow.available
        ? "<span>Ready to run.</span>"
        : `<span class="missing">${escapeHtml(workflow.executable)} is not available on PATH.</span>`
    }
  `;
}

async function loadWorkflows() {
  const response = await fetch("/api/workflows");
  const payload = await response.json();
  workflows = payload.workflows;
  workflowSelect.innerHTML = workflows
    .map((workflow) => `<option value="${escapeHtml(workflow.id)}">${escapeHtml(workflow.name)}</option>`)
    .join("");
  renderWorkflow();
}

function renderFiles(files) {
  filesBox.innerHTML = "";
  for (const file of files || []) {
    const link = document.createElement("a");
    link.href = `/${file.replaceAll("\\", "/")}`;
    link.target = "_blank";
    link.rel = "noreferrer";
    link.textContent = file;
    filesBox.appendChild(link);
  }
}

function citationLink(citation) {
  const title = escapeHtml(citation.title || citation.source || "Citation");
  const meta = [citation.source, citation.year, citation.pmid ? `PMID ${citation.pmid}` : null]
    .filter(Boolean)
    .map(escapeHtml)
    .join(" · ");
  const label = meta ? `${title}<span>${meta}</span>` : title;
  if (citation.url) {
    return `<a href="${escapeHtml(citation.url)}" target="_blank" rel="noreferrer">${label}</a>`;
  }
  return `<span>${label}</span>`;
}

function renderInhibitors(payload) {
  inhibitorCount.textContent = `${payload.count} found`;
  inhibitorStatus.textContent = payload.disclaimer;
  targetCard.classList.remove("hidden");
  targetCard.innerHTML = `
    <strong>${escapeHtml(payload.target.name)}</strong>
    <span>${escapeHtml(payload.target.organism)} · ${escapeHtml(payload.target.type || "target")}</span>
    <a href="${escapeHtml(payload.target.url)}" target="_blank" rel="noreferrer">${escapeHtml(payload.target.chembl_id)}</a>
    <small>Matched symbols: ${escapeHtml((payload.target.symbols || []).join(", ") || payload.gene)}</small>
  `;

  if (!payload.results.length) {
    inhibitorList.innerHTML = `
      <article class="empty-state">
        <h3>No inhibitor mechanisms found</h3>
        <p>No exact matching inhibitor evidence was returned. If you supplied a disease/tissue context, this means no open-access citation explicitly matched that exact context for the target and inhibitor.</p>
      </article>
    `;
    return;
  }

  inhibitorList.innerHTML = payload.results
    .map((item, index) => {
      const mechanisms = (item.mechanisms || [])
        .map((mechanism) => `
          <li>
            <strong>${escapeHtml(mechanism.action_type || "Mechanism")}</strong>
            ${escapeHtml(mechanism.mechanism_of_action || "No mechanism text available")}
            ${mechanism.comment ? `<span>${escapeHtml(mechanism.comment)}</span>` : ""}
          </li>
        `)
        .join("");
      const citations = (item.citations || [])
        .map((citation) => `<li>${citationLink(citation)}<em>${escapeHtml(citation.evidence_source || "source")}</em></li>`)
        .join("");
      return `
        <article class="inhibitor-card">
          <header>
            <span class="rank">${index + 1}</span>
            <div>
              <h3>${escapeHtml(item.name)}</h3>
              <a href="${escapeHtml(item.chembl_url)}" target="_blank" rel="noreferrer">${escapeHtml(item.chembl_id)}</a>
            </div>
            <div class="phase">${item.max_phase ? `Phase ${escapeHtml(item.max_phase)}` : "Phase unknown"}</div>
          </header>
          <dl>
            <div><dt>First approval</dt><dd>${escapeHtml(item.first_approval || "Not listed")}</dd></div>
            <div><dt>Type</dt><dd>${escapeHtml(item.molecule_type || "Small molecule / not listed")}</dd></div>
          </dl>
          <h4>Mechanism evidence</h4>
          <ul class="mechanisms">${mechanisms}</ul>
          <h4>Citations</h4>
          ${item.context_match_required ? `<p class="context-badge">${escapeHtml(item.context_citation_count)} exact-context citation(s)</p>` : ""}
          <ul class="citations">${citations || "<li>No citation links returned for this mechanism record.</li>"}</ul>
        </article>
      `;
    })
    .join("");
}

tabs.forEach((tab) => {
  tab.addEventListener("click", () => {
    tabs.forEach((item) => item.classList.remove("active"));
    panels.forEach((panel) => panel.classList.remove("active"));
    tab.classList.add("active");
    document.querySelector(`#${tab.dataset.panel}`).classList.add("active");
  });
});

inhibitorForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  inhibitorButton.disabled = true;
  inhibitorButton.textContent = "Searching...";
  inhibitorStatus.textContent = "Searching ChEMBL, PubMed, and Europe PMC open records...";
  inhibitorCount.textContent = "Working";
  targetCard.classList.add("hidden");
  inhibitorList.innerHTML = "";

  const gene = document.querySelector("#gene-query").value.trim();
  const context = document.querySelector("#context-query").value.trim();
  const exclude = document.querySelector("#exclude-query").value.trim();
  const requireContext = document.querySelector("#require-context").checked ? "true" : "false";
  const limit = document.querySelector("#inhibitor-limit").value || "50";
  const params = new URLSearchParams({ gene, limit, context, exclude, require_context: requireContext });

  try {
    const response = await fetch(`/api/inhibitors?${params.toString()}`);
    const contentType = response.headers.get("content-type") || "";
    const payload = contentType.includes("application/json")
      ? await response.json()
      : { error: `Server returned ${response.status} ${response.statusText || "non-JSON response"}. Restart python app.py, then refresh this page.` };
    if (!response.ok || !payload.ok) {
      throw payload;
    }
    renderInhibitors(payload);
  } catch (error) {
    inhibitorCount.textContent = "Error";
    inhibitorStatus.textContent = error.error || error.message || "Search failed.";
    inhibitorList.innerHTML = "";
  } finally {
    inhibitorButton.disabled = false;
    inhibitorButton.textContent = "Find inhibitors";
  }
});

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  button.disabled = true;
  button.textContent = "Running...";
  statusBox.textContent = "Workflow started.";
  stdoutBox.textContent = "";
  filesBox.innerHTML = "";

  const parameters = {
    threads: document.querySelector("#threads").value,
    database: document.querySelector("#database").value,
  };

  const body = new FormData(form);
  body.append("parameters", JSON.stringify(parameters));

  try {
    const response = await fetch("/api/run", {
      method: "POST",
      body,
    });
    const payload = await response.json();
    if (!response.ok || !payload.ok) {
      throw payload;
    }
    statusBox.textContent = `Completed ${payload.workflow} in run ${payload.run_id}.`;
    stdoutBox.textContent = payload.stdout || JSON.stringify(payload.summary, null, 2);
    renderFiles(payload.files);
  } catch (error) {
    statusBox.textContent = `Run failed: ${error.error || error.message || "Unknown error"}`;
    stdoutBox.textContent = JSON.stringify(error, null, 2);
    renderFiles(error.files);
  } finally {
    button.disabled = false;
    button.textContent = "Run workflow";
  }
});

workflowSelect.addEventListener("change", renderWorkflow);
loadWorkflows();
