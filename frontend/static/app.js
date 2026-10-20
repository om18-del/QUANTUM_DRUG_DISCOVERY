const API_BASE = (() => {
  if (window.__QD_API_BASE__) return window.__QD_API_BASE__;
  return "/api";
})();

// The stylesheet is injected into index.html by the server (app/main.py) with a
// version query string, and served no-cache. Do not add a second <link> here: a
// hard-coded version such as "?v=7" diverges from the server's, and because a
// later-appended sheet wins the cascade it silently reintroduces the stale-CSS
// problem that the server-side injection exists to prevent.

const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));

function setStatus(el, kind, message) {
  if (!el) return;
  el.className = "status-line";
  el.innerHTML = "";
  if (message == null || message === "") {
    const span = document.createElement("span");
    span.className = "muted";
    span.textContent = "—";
    el.appendChild(span);
    return;
  }
  const span = document.createElement("span");
  span.className = String(kind || "").trim() || "muted";
  span.textContent = String(message);
  el.appendChild(span);
}

function setBadge(kind) {
  const el = $("#app-status");
  if (!el) return;
  const k = String(kind || "").trim();
  let label = "connecting";
  if (k === "ok") label = "api ok";
  else if (k === "warn") label = "api warning";
  else if (k === "err") label = "api error";
  else if (k === "placeholder") label = "placeholder";
  else if (k === "offline") label = "offline";
  el.className = "status-badge status-" + (k || "idle");
  el.textContent = label;
}

async function apiFetch(path, init = {}) {
  const url = API_BASE + path;
  let response;
  try {
    response = await fetch(url, init);
  } catch (err) {
    console.warn("apiFetch network error", url, err);
    const error = { code: "network_error", message: err && err.message ? err.message : "Network request failed" };
    throw { status: 0, error };
  }
  let body = null;
  try {
    body = await response.json();
  } catch (err) {
    console.warn("apiFetch non-json response", url, response.status, err);
    body = { error: { code: "bad_response", message: "Server returned a non-JSON response" } };
  }
  if (!response.ok) {
    const error = (body && body.error) ? body.error : { code: "http_error", message: "Request failed" };
    throw { status: response.status, error };
  }
  if (body && body.error) {
    console.warn("apiFetch upstream error", url, body.error);
  }
  return { ok: true, status: response.status, body };
}

function escapeHtml(s) {
  if (s == null) return "";
  if (typeof s === "string") {
    return s
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }
  return String(s);
}

function errorNote(error) {
  const code = error && error.code ? error.code : "unknown";
  const message = error && error.message ? error.message : String(error);
  return `<div class="error-note"><strong>${escapeHtml(code)}</strong>: ${escapeHtml(message)}</div>`;
}

function placeholderNote(note) {
  if (!note) return "";
  return `<div class="placeholder-note">${escapeHtml(note)}</div>`;
}

/**
 * Collapsed one-click detail for explanations that are too long to lead with.
 * File paths and Colab instructions belong here (and in the README), not in the
 * panel where they push the buttons around.
 */
function detailNote(summary, detail) {
  if (!detail) return "";
  return `<details class="detail-note"><summary>${escapeHtml(summary)}</summary><div>${escapeHtml(detail)}</div></details>`;
}

/**
 * One short plain sentence for an unavailable panel, with the long technical
 * note (file paths, Colab instructions) collapsed underneath it.
 */
function unavailableBlock(shortText, longNote) {
  return `<div class="placeholder-note">${escapeHtml(shortText)}</div>` + detailNote("Why this is unavailable", longNote);
}

/**
 * Short state for a status line. The full note or error is rendered inside the
 * result card, so repeating it here would show the same message twice.
 */
function shortState(body) {
  if (!body) return "failed";
  if (body.status === "computed") return "computed";
  if (body.status === "placeholder") return "not available";
  return String(body.status || "unknown");
}

function shortErrorState(error) {
  const code = error && error.code ? String(error.code) : "error";
  if (code === "compound_not_found") return "no match";
  if (code === "network_error") return "no connection";
  if (code === "not_found") return "not found";
  if (code === "hamiltonian_unavailable") return "not available";
  return "failed";
}

/* Short, plain "not available" wording. Every panel uses this same sentence
   style so an unavailable panel looks the same as the others. The resource
   estimator needs no entry here: it counts from the formula and never depends
   on an artifact being present. */
const NOT_AVAILABLE = {
  curve: "The H2 energy curve has not been generated yet.",
  simulate: "H2 simulation is not available yet.",
};

/** Numbers arrive as JSON numbers or numeric strings; never print "null". */
function num(value, digits = 2) {
  if (value == null || value === "") return "—";
  const n = Number(value);
  if (!isFinite(n)) return escapeHtml(value);
  return n.toFixed(digits);
}

function formatEnergy(ha) {
  if (ha == null || ha === "") return "—";
  return ha.toFixed(4) + " Ha";
}

function formatMha(mha) {
  if (mha == null) return "—";
  return mha.toFixed(2) + " mHa";
}

function jsonBlock(obj, extraClasses = "") {
  return `<details><summary>JSON</summary><pre class="${extraClasses}">${escapeHtml(JSON.stringify(obj, null, 2))}</pre></details>`;
}

function showPanel(name) {
  $$(".tab").forEach((t) => t.classList.toggle("active", t.dataset.tab === name));
  $$(".panel").forEach((p) => p.hidden = p.dataset.panel !== name);
}

function initTabs() {
  const hash = location.hash.replace("#", "") || "quantum";
  if ($$(".tab[data-tab='" + hash + "']").length) {
    showPanel(hash);
  } else {
    showPanel("quantum");
  }
  $$(".tab").forEach((t) =>
    t.addEventListener("click", () => {
      showPanel(t.dataset.tab);
      history.replaceState(null, "", "#" + t.dataset.tab);
    })
  );
  window.addEventListener("hashchange", () => {
    const target = location.hash.replace("#", "");
    if ($$(".tab[data-tab='" + target + "']").length) {
      showPanel(target);
    }
  });
}

function clearStatus(selector) {
  const el = $(selector);
  if (!el) return;
  el.className = "status-line";
  el.innerHTML = "";
}

function clearPanelStatuses() {
  [
    "#h2-curve-status",
    "#simulate-status",
    "#resource-status",
    "#search-status",
    "#props-status",
    "#sim-status",
    "#protein-status",
    "#complexes-status",
    "#viewer-status",
    "#protein-viewer-status",
    "#complex-viewer-status",
  ].forEach(clearStatus);
}

/* ---------------------------------------------------------------- 3Dmol --- */

/* Jmol/CPK element colours. 3Dmol ships this palette under the name "Jmol";
   naming it explicitly means a style change can never silently change it. */
const ELEMENT_SCHEME = "Jmol";

/* Start tilted so the depth reads immediately instead of a face-on view that
   could be mistaken for a flat drawing. */
const TILT_Y_DEG = 27;
const TILT_X_DEG = 25;

/* Gentle auto-rotate, not a spinning top. */
const SPIN_SPEED = 0.3;

/* Labelling every atom of a protein would be unreadable, so the label toggle
   is capped and, above the cap, says how many atoms it left out. */
const ATOM_LABEL_LIMIT = 60;

/* Cartoon colours for chains, cycled over however many chains are present. */
const CHAIN_COLORS = ["#6ea8fe", "#f2c94a", "#6fcf97", "#ff8f6b", "#c4a2f5", "#7c6bfe", "#ff6b6b", "#5ad1e6"];
const WATER_NAMES = ["HOH", "WAT", "DOD", "H2O", "SOL"];

function getViewerLib() {
  // 3Dmol's UMD bundle exposes the global as window["3Dmol"].
  return window["3Dmol"] || window.$3Dmol || window.$_3Dmol || null;
}

function createViewer(el) {
  const lib = getViewerLib();
  if (!lib || typeof lib.createViewer !== "function") return null;
  el.innerHTML = "";
  // 3Dmol clears its canvas to an opaque colour, so the gradient cannot live
  // behind it; the darkest stop is used here and the lift is painted over the
  // canvas by .viewer-wrap::after in the stylesheet.
  const viewer = lib.createViewer(el, { backgroundColor: "#05070d", antialias: true });
  try {
    if (typeof viewer.setCameraType === "function") viewer.setCameraType("perspective");
  } catch (e) { /* older builds already default to perspective */ }
  return viewer;
}

/** Zoom to fit, then tilt. Order matters: zoomTo resets the orientation. */
function frameAndTilt(viewer) {
  try { viewer.zoomTo(); } catch (e) { /* ignore */ }
  try {
    if (typeof viewer.rotate === "function") {
      viewer.rotate((TILT_Y_DEG * Math.PI) / 180, "y");
      viewer.rotate((TILT_X_DEG * Math.PI) / 180, "x");
    }
  } catch (e) { /* a build without rotate still renders, just face-on */ }
}

function spinViewer(viewer, on) {
  try { viewer.spin(on ? SPIN_SPEED : false); } catch (e) { /* ignore */ }
}

/* Shared presets, so "ball-and-stick" means the same thing in all three
   viewers: thin sticks, spheres at 0.28 of the van der Waals radius. */
function ballAndStickSpec() {
  return {
    stick: { radius: 0.12, colorscheme: ELEMENT_SCHEME },
    sphere: { scale: 0.28, colorscheme: ELEMENT_SCHEME },
  };
}

function styleSpec(style) {
  if (style === "stick") return { stick: { radius: 0.12, colorscheme: ELEMENT_SCHEME } };
  if (style === "spacefill") return { sphere: { scale: 1.0, colorscheme: ELEMENT_SCHEME } };
  return ballAndStickSpec();
}

/* Hydrogens drawn at the same scale as the heavy atoms they hang off swamp the
   structure, so they get their own, smaller, style. */
function addHydrogenStyle(viewer) {
  try {
    viewer.addStyle({ elem: "H" }, {
      sphere: { scale: 0.16, colorscheme: ELEMENT_SCHEME },
      stick: { radius: 0.09, colorscheme: ELEMENT_SCHEME },
    });
  } catch (e) { /* ignore */ }
}

function translucentSurface(viewer) {
  const lib = getViewerLib();
  const type = lib && lib.SurfaceType ? lib.SurfaceType.VDW : 1;
  try { viewer.addSurface(type, { opacity: 0.42, color: "0x6ea8fe" }, {}); } catch (e) { /* optional */ }
}

/** Replace the tracked element labels, leaving any caller's labels alone. */
function setElementLabels(viewer, atoms, store) {
  (store || []).forEach((l) => { try { viewer.removeLabel(l); } catch (e) { /* ignore */ } });
  if (store) store.length = 0;
  if (!atoms || !store) return;
  const shown = atoms.slice(0, ATOM_LABEL_LIMIT);
  shown.forEach((a) => {
    const text = a.elem || a.symbol || "";
    if (!text) return;
    store.push(viewer.addLabel(text, {
      position: a,
      fontColor: "#e6eaf4",
      backgroundColor: "#0b0f1c",
      backgroundOpacity: 0.5,
      fontSize: 10,
      showBackground: false,
    }));
  });
}

function allAtoms(viewer) {
  try {
    const model = viewer.getModel();
    return model && typeof model.selectedAtoms === "function" ? model.selectedAtoms({}) : [];
  } catch (e) {
    return [];
  }
}

/** Non-water hetero residues — the bound ligand(s), not the crystallisation solvent. */
function ligandSelection(viewer) {
  try {
    const atoms = viewer.getModel().selectedAtoms({ hetflag: true });
    const names = [];
    atoms.forEach((a) => {
      const n = String(a.resn || "").toUpperCase();
      if (WATER_NAMES.indexOf(n) === -1 && names.indexOf(n) === -1) names.push(n);
    });
    return names.length ? { hetflag: true, resn: names } : null;
  } catch (e) {
    return null;
  }
}

/** A5: bound ligands are always drawn as ball-and-stick, whatever the protein is. */
function addLigandSticks(viewer) {
  const spec = ligandSelection(viewer);
  if (!spec) return;
  try { viewer.addStyle(spec, ballAndStickSpec()); } catch (e) { /* ignore */ }
}

function chainsOf(viewer) {
  const seen = [];
  allAtoms(viewer).forEach((a) => {
    const c = a.chain == null ? "" : String(a.chain);
    if (seen.indexOf(c) === -1) seen.push(c);
  });
  return seen.sort();
}

/** A5: colour the cartoon by chain so two chains of one protein read apart. */
function colourByChain(viewer) {
  const chains = chainsOf(viewer);
  if (chains.length < 2) {
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    return;
  }
  viewer.setStyle({}, { cartoon: { color: CHAIN_COLORS[0] } });
  chains.forEach((chain, i) => {
    viewer.setStyle({ chain: chain }, { cartoon: { color: CHAIN_COLORS[i % CHAIN_COLORS.length] } });
  });
}

/** Protein-or-ligand styling shared by the Protein Viewer and the Studio. */
function applyProteinStyle(viewer, opts) {
  const options = opts || {};
  const style = options.style || "cartoon";
  try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
  if (style === "cartoon") colourByChain(viewer);
  else viewer.setStyle({}, styleSpec(style));
  addLigandSticks(viewer);
  setElementLabels(viewer, options.labels ? allAtoms(viewer) : null, options.labelStore);
}

function toText(data) {
  if (data == null) return "";
  if (typeof data === "string") return data;
  if (data instanceof ArrayBuffer) return new TextDecoder().decode(data);
  if (typeof data.length === "number") {
    try {
      return new TextDecoder().decode(new Uint8Array(data));
    } catch (e) {
      return String(data);
    }
  }
  return String(data);
}

/** Centre of mass of a 3Dmol selection, in model coordinates. */
function centroid(viewer, spec) {
  try {
    const model = viewer.getModel();
    const atoms = model && typeof model.selectedAtoms === "function" ? model.selectedAtoms(spec) : [];
    if (!atoms || !atoms.length) return null;
    let x = 0, y = 0, z = 0;
    atoms.forEach((a) => { x += a.x; y += a.y; z += a.z; });
    return { x: x / atoms.length, y: y / atoms.length, z: z / atoms.length };
  } catch (e) {
    return null;
  }
}

function residueSpec(residue) {
  const spec = { resi: String(residue.residue_number) };
  if (residue.chain) spec.chain = residue.chain;
  return spec;
}

const HIGHLIGHT = { stick: { radius: 0.28, color: "#ffd166" }, sphere: { radius: 0.42, color: "#ffd166" } };
const CATALYTIC_HIGHLIGHT = { stick: { radius: 0.3, color: "#ff6b6b" }, sphere: { radius: 0.45, color: "#ff6b6b" } };

window.addEventListener("DOMContentLoaded", () => {
  $("#api-base-label").textContent = "api: " + API_BASE;
  $("#foot-meta").textContent = `serving from ${API_BASE} · every figure is computed here or labelled as literature`;
  initTabs();
  setBadge("idle");
  clearPanelStatuses();

  const onClick = (sel, fn) => { const el = $(sel); if (el) el.addEventListener("click", fn); };
  const onChange = (sel, fn) => { const el = $(sel); if (el) el.addEventListener("change", fn); };

  onClick("#refresh-curve", loadCurve);
  onClick("#refresh-hadamard", loadHadamard);
  onClick("#simulate-run", runSimulate);
  onClick("#resource-run", runResources);
  onClick("#search-run", runSearch);
  onClick("#props-run", runProperties);
  onClick("#sim-run", runMoleculeSimulation);
  onClick("#protein-fetch", fetchProtein);
  onClick("#complexes-list", onComplexClick);
  onClick("#explore-run", () => runExplore($("#explore-pdb").value));
  onClick("#explore-candidates", onExploreCandidateClick);
  onChange("#complex-category", () => {
    complexesState.category = $("#complex-category").value;
    renderComplexList();
  });
  const search = $("#complex-search");
  if (search) {
    search.addEventListener("input", () => {
      complexesState.query = search.value;
      renderComplexList();
    });
  }
  const exploreInput = $("#explore-pdb");
  if (exploreInput) {
    exploreInput.addEventListener("keydown", (event) => {
      if (event.key === "Enter") {
        event.preventDefault();
        runExplore(exploreInput.value);
      }
    });
  }

  onClick("#mol-reset", () => resetMolView());
  onClick("#mol-2d", show2DDepiction);
  onChange("#mol-style", () => redrawMolecule());
  onChange("#mol-labels", () => redrawMolecule());
  onChange("#mol-rotate", () => toggleSpin(molView, $("#mol-rotate").checked));
  onChange("#mol-conformer", () => showConformer(parseInt($("#mol-conformer").value, 10) || 0));

  onChange("#protein-style", () => redrawProtein());
  onChange("#protein-labels", () => redrawProtein());
  onChange("#protein-rotate", () => toggleSpin(proteinView, $("#protein-rotate").checked));

  onClick("#complex-play", playMechanism);
  onClick("#complex-reset", () => resetComplexView());
  onChange("#complex-style", () => redrawComplex());
  onChange("#complex-labels", () => redrawComplex());
  onChange("#complex-rotate", () => toggleSpin(complexView, $("#complex-rotate").checked));

  (async () => {
    try {
      await Promise.allSettled([
        loadCurve(),
        loadHadamard(),
        loadComplexes(),
      ]);
      setBadge("ok");
    } catch (e) {
      console.warn("startup failed", e);
      setBadge("err");
    }
  })();
});

/* --------------------------------------------------------------- Quantum --- */

/* The last simulate run's bond-length disclosure, so the curve can mark where the
   point came from. Null until a run has succeeded. */
let lastSimulateMarker = null;

/* The last successfully loaded curve payload, so a simulate run can re-mark the
   existing chart without another network round-trip. */
let lastCurveBody = null;

async function loadCurve() {
  const canvas = $("#h2-curve");
  const status = $("#h2-curve-status");
  const meta = $("#curve-meta");
  const empty = $("#h2-curve-empty");
  const detail = $("#h2-curve-detail");

  // Same short wording the other panels use, with the long file-path/Colab
  // explanation collapsed underneath instead of pushed into the layout.
  const showUnavailable = (short, longText) => {
    if (empty) empty.textContent = short;
    if (detail) detail.innerHTML = detailNote("Why this is unavailable", longText);
    setStatus(status, "placeholder", "not available");
    meta.textContent = "";
    lastCurveBody = null;
  };

  if (!canvas || !canvas.getContext) {
    setStatus(status, "err", "canvas unavailable");
    return;
  }
  const ctx = canvas.getContext("2d");
  canvas.width = Math.max(2, canvas.clientWidth || 300);
  canvas.height = Math.max(2, canvas.clientHeight || 220);
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  try {
    const { body } = await apiFetch("/quantum/h2-curve");
    if (body.status !== "computed") {
      showUnavailable(NOT_AVAILABLE.curve, body.note);
      return;
    }
    if (empty) empty.textContent = "";
    if (detail) detail.innerHTML = "";
    setStatus(status, "ok", "curve computed");
    meta.textContent = `${body.molecule || "—"} · ${body.basis || "—"} · ${body.bond_lengths_angstrom ? body.bond_lengths_angstrom.length : 0} points`;
    lastCurveBody = body;
    drawCurve(ctx, canvas, body);
  } catch (e) {
    const err = e && e.error ? e.error : e;
    if (err && err.code === "hamiltonian_unavailable") {
      // Expected state, not a fault: report it as "not available" like every
      // other panel, and keep the long explanation one click away.
      showUnavailable(NOT_AVAILABLE.curve, err.message);
      return;
    }
    if (empty) empty.textContent = "The curve could not be loaded.";
    if (detail) detail.innerHTML = "";
    setStatus(status, "err", shortErrorState(err));
    meta.textContent = "";
  }
}

function drawCurve(ctx, canvas, body) {
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  const cssW = Math.max(1, canvas.clientWidth || 300);
  const cssH = Math.max(1, canvas.clientHeight || 220);
  canvas.width = Math.floor(cssW * dpr);
  canvas.height = Math.floor(cssH * dpr);
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.scale(dpr, dpr);
  const pad = { top: 14, right: 14, bottom: 40, left: 44 };
  const plotW = Math.max(10, cssW - pad.left - pad.right);
  const plotH = Math.max(10, cssH - pad.top - pad.bottom);

  ctx.fillStyle = "#0b0f1c";
  ctx.fillRect(0, 0, cssW, cssH);

  ctx.textBaseline = "alphabetic";

  const xs = body.bond_lengths_angstrom;
  const exact = body.exact_energies_ha;
  const vqe = body.vqe_energies_ha;

  const minE = Math.min(...exact);
  const maxE = Math.max(...exact);
  const range = Math.max(1e-9, maxE - minE);
  const xScale = (i) => pad.left + (i / Math.max(1, xs.length - 1)) * plotW;
  const yScale = (v) => pad.top + plotH - ((v - minE) / range) * plotH;

  ctx.strokeStyle = "#8892b5";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, pad.top + plotH);
  ctx.lineTo(pad.left + plotW, pad.top + plotH);
  ctx.stroke();

  ctx.fillStyle = "#8892b5";
  ctx.font = "10px sans-serif";
  ctx.textAlign = "right";
  for (let i = 0; i <= 4; i++) {
    const v = minE + (range * i) / 4;
    const y = yScale(v);
    ctx.fillText(v.toFixed(3), pad.left - 6, y + 3);
    ctx.strokeStyle = "#1c2337";
    ctx.beginPath();
    ctx.moveTo(pad.left, y);
    ctx.lineTo(pad.left + plotW, y);
    ctx.stroke();
  }
  // X ticks: six labelled values, not one label per stored point. All 58 point
  // labels used to overlap into an unreadable smear. Every value below sits on
  // the stored grid (0.20 .. 3.00 in 0.05 steps, plus 0.74), and a value that is
  // not on it is skipped rather than drawn at the wrong place.
  const xTickValues = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
  const gridStep = xs.length > 1 ? Math.abs(xs[1] - xs[0]) : 0.05;
  const xTolerance = Math.max(1e-9, gridStep / 2);
  ctx.strokeStyle = "#8892b5";
  ctx.fillStyle = "#8892b5";
  ctx.textAlign = "center";
  for (const value of xTickValues) {
    let nearestIndex = -1;
    let nearestDistance = Infinity;
    for (let i = 0; i < xs.length; i++) {
      const distance = Math.abs(xs[i] - value);
      if (distance < nearestDistance) {
        nearestDistance = distance;
        nearestIndex = i;
      }
    }
    if (nearestIndex < 0 || nearestDistance > xTolerance) continue;
    const x = xScale(nearestIndex);
    ctx.beginPath();
    ctx.moveTo(x, pad.top + plotH);
    ctx.lineTo(x, pad.top + plotH + 4);
    ctx.stroke();
    ctx.fillText(value.toFixed(1), x, pad.top + plotH + 17);
  }
  // The unit is stated once, under the ticks, instead of once per point.
  ctx.fillText("bond length (Å)", pad.left + plotW / 2, pad.top + plotH + 33);

  if (exact) {
    ctx.strokeStyle = "#6ea8fe";
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < exact.length; i++) {
      const x = xScale(i);
      const y = yScale(exact[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  if (vqe) {
    ctx.strokeStyle = "#f2c94a";
    ctx.setLineDash([4, 4]);
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < vqe.length; i++) {
      if (vqe[i] == null) continue;
      const x = xScale(i);
      const y = yScale(vqe[i]);
      if (i === 0 || vqe[i - 1] == null) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  ctx.fillStyle = "#e6eaf4";
  ctx.font = "11px sans-serif";
  ctx.textAlign = "left";
  const legX = pad.left + 8;
  const legY = pad.top + 14;
  ctx.fillStyle = "#6ea8fe";
  ctx.fillRect(legX, legY - 5, 12, 2);
  ctx.fillStyle = "#e6eaf4";
  ctx.fillText("exact", legX + 16, legY);
  ctx.fillStyle = "#f2c94a";
  ctx.fillRect(legX + 80, legY - 5, 12, 2);
  ctx.fillStyle = "#e6eaf4";
  ctx.fillText("VQE", legX + 96, legY);

  // Mark the bond length the last simulate run used: a vertical guide, a dot on
  // the exact curve, and the snap disclosure ("requested X, used Y") when the
  // request was off-grid. Nothing is drawn until a run has succeeded.
  if (lastSimulateMarker && Array.isArray(xs) && xs.length && Array.isArray(exact) && exact.length) {
    const used = Number(lastSimulateMarker.used);
    if (Number.isFinite(used)) {
      let nearestIndex = 0;
      let nearestDist = Infinity;
      for (let i = 0; i < xs.length; i++) {
        const d = Math.abs(xs[i] - used);
        if (d < nearestDist) { nearestDist = d; nearestIndex = i; }
      }
      const mx = xScale(nearestIndex);
      const my = yScale(exact[nearestIndex]);
      ctx.strokeStyle = "#ff8f6b";
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(mx, pad.top);
      ctx.lineTo(mx, pad.top + plotH);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = "#ff8f6b";
      ctx.beginPath();
      ctx.arc(mx, my, 3.5, 0, Math.PI * 2);
      ctx.fill();

      const requested = Number(lastSimulateMarker.requested);
      const label = lastSimulateMarker.snapped && Number.isFinite(requested)
        ? `used ${used.toFixed(3)} Å (requested ${requested.toFixed(3)})`
        : `used ${used.toFixed(3)} Å`;
      ctx.font = "10px sans-serif";
      ctx.textAlign = "left";
      const labelWidth = Math.min(ctx.measureText(label).width + 8, cssW - mx - 6);
      ctx.fillStyle = "rgba(11, 15, 28, 0.86)";
      ctx.fillRect(mx + 4, pad.top + plotH - 15, Math.max(0, labelWidth), 14);
      ctx.fillStyle = "#ff8f6b";
      ctx.fillText(label, mx + 7, pad.top + plotH - 4);
    }
  }
}

/** Re-render the cached curve payload, e.g. to add a new simulate marker. */
function redrawCurveWithMarker() {
  const canvas = $("#h2-curve");
  if (!canvas || !canvas.getContext || !lastCurveBody) return;
  drawCurve(canvas.getContext("2d"), canvas, lastCurveBody);
}

async function loadHadamard() {
  const el = $("#hadamard-result");
  el.hidden = false;
  try {
    const { body } = await apiFetch("/quantum/hadamard-demo");
    if (body.status !== "computed") {
      el.innerHTML = placeholderNote(body.status ? `demo status: ${body.status}` : "unknown demo status") + jsonBlock(body);
      return;
    }
    el.innerHTML = `
      <div class="computed-note">computed</div>
      <dl>
        <dt>circuit</dt><dd>${escapeHtml(body.circuit)}</dd>
        <dt>initial state</dt><dd>[${body.initial_state}]</dd>
        <dt>final state</dt><dd>[${body.final_state ? body.final_state.map(v => Number(v).toFixed(4)).join(", ") : "—"}]</dd>
        <dt>P(0)</dt><dd>${body.probabilities && body.probabilities["0"] != null ? Number(body.probabilities["0"]).toFixed(4) : "—"}</dd>
        <dt>P(1)</dt><dd>${body.probabilities && body.probabilities["1"] != null ? Number(body.probabilities["1"]).toFixed(4) : "—"}</dd>
      </dl>
      ${jsonBlock(body)}
    `;
  } catch (e) {
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

async function runSimulate() {
  const bond = $("#sim-bond").value;
  const molecule = $("#sim-molecule").value;
  // Default run mode is the exact statevector (shots: null); the dropdown's value
  // is either "statevector" or a shot count. There is no free-text shots box.
  const runMode = $("#sim-run-mode") ? $("#sim-run-mode").value : "statevector";
  const shots = runMode === "statevector" ? null : parseInt(runMode, 10);
  // The noise box is disabled in the form and labelled "not applied yet": the
  // backend accepts and echoes the field but no noise model exists, so the value
  // is inert. It is still sent so the request shape stays the documented one.
  const noise = $("#sim-noise").value;
  const status = $("#simulate-status");
  const el = $("#simulate-result");
  el.hidden = false;
  setStatus(status, "placeholder", "running…");

  try {
    const { body } = await apiFetch("/quantum/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ bond_length: parseFloat(bond), shots, noise: parseFloat(noise), molecule }),
    });

    if (body.status === "placeholder") {
      setStatus(status, "placeholder", "not available");
      el.innerHTML = unavailableBlock(NOT_AVAILABLE.simulate, body.note) + jsonBlock(body);
      return;
    }
    if (body.status === "failed") {
      setStatus(status, "err", "failed");
      el.innerHTML = `<div class="error-note">simulation failed</div>` + jsonBlock(body);
      return;
    }

    // The panel marks where this run's energy came from, and says so when the
    // requested bond length was off the stored grid.
    const used = Number(body.used_bond_length ?? body.bond_length);
    const requested = Number(body.requested_bond_length ?? body.bond_length);
    const snapped = Boolean(body.snapped) || (Number.isFinite(requested) && Number.isFinite(used) && Math.abs(requested - used) > 1e-9);
    lastSimulateMarker = Number.isFinite(used)
      ? { used, requested, snapped }
      : null;
    redrawCurveWithMarker();

    const isShot = body.shots != null;
    const stdHa = body.vqe_energy_std_ha;
    const budget = `${body.n_restarts ?? body.restarts ?? "—"} restarts × ${body.total_energy_evaluations ?? body.iterations ?? "—"} evaluations`;
    const snapRow = snapped
      ? `<dt>bond-length snap</dt><dd>requested ${requested.toFixed(3)} Å, used ${used.toFixed(3)} Å</dd>`
      : "";
    const uncertaintyRow = isShot && stdHa != null
      ? `<dt>VQE uncertainty</dt><dd>± ${(Number(stdHa) * 1000).toFixed(2)} mHa (${escapeHtml(String(body.vqe_energy_uncertainty_method || "").replace(/_/g, " "))})</dd>`
      : "";
    // The 1.6 mHa chemical-accuracy badge is a noiseless claim only; a shot run gets
    // the sampling-error note instead, never a pass/fail badge.
    const accuracyBadge = isShot
      ? `<div class="shot-note">shot noise: sampling error, not a hardware result</div>`
      : (body.within_chemical_accuracy
          ? `<div class="accuracy-badge ok">within 1.6 mHa of exact</div>`
          : `<div class="accuracy-badge warn">outside the 1.6 mHa target</div>`);

    setStatus(status, "ok", `${body.molecule} · ${used.toFixed(3)} Å`);
    el.innerHTML = `
      <div class="computed-note">status: ${escapeHtml(body.status)}</div>
      <dl>
        <dt>molecule</dt><dd>${escapeHtml(body.molecule)}</dd>
        <dt>bond length used</dt><dd>${used.toFixed(3)} Å</dd>
        ${snapRow}
        <dt>HF energy</dt><dd>${formatEnergy(body.hf_energy_ha)}</dd>
        <dt>exact energy</dt><dd>${formatEnergy(body.exact_energy_ha)}</dd>
        <dt>VQE energy</dt><dd>${formatEnergy(body.vqe_energy_ha)}</dd>
        ${uncertaintyRow}
        <dt>error vs exact</dt><dd>${formatMha(body.error_mha)}</dd>
        <dt>run mode</dt><dd>${isShot ? `${body.shots} shots` : "exact statevector"}</dd>
        <dt>optimisation budget</dt><dd>${budget}</dd>
        <dt>estimator</dt><dd>${escapeHtml(body.estimator_class_name || "—")}</dd>
      </dl>
      ${accuracyBadge}
      ${jsonBlock(body)}
    `;
  } catch (e) {
    const err = e && e.error ? e.error : e;
    if (err && err.code === "hamiltonian_unavailable") {
      // Same expected state the curve card reports, so it gets the same badge
      // style and the same short sentence: a missing artifact makes the energy
      // unavailable, it is not a fault.
      setStatus(status, "placeholder", "not available");
      el.innerHTML = unavailableBlock(NOT_AVAILABLE.simulate, err.message);
      return;
    }
    setStatus(status, "err", shortErrorState(err));
    el.innerHTML = errorNote(err);
  }
}

/* ------------------------------------------------------- Resource estimator --- */

/* The resource estimator answers "what would this molecule need?" from a formula
   alone, and every figure it returns is tagged in the backend's `tags` dict as
   exact (a closed-form count) or estimated (a judgement or a rough limit). The
   table below prints those tags rather than inventing its own wording. */

const RESOURCE_FIXED_NOTE = "Computed by formula for the STO-3G basis. No simulation was run for this molecule.";

function tagChip(tag) {
  if (!tag) return "";
  return `<span class="chip ${escapeHtml(tag)}">${escapeHtml(tag)}</span>`;
}

function quantityRows(pairs) {
  return pairs
    .map(([label, value, tag]) => `
      <tr>
        <td>${escapeHtml(label)}</td>
        <td class="numeric">${escapeHtml(value == null ? "—" : String(value))}</td>
        <td>${tagChip(tag)}</td>
      </tr>
    `)
    .join("");
}

function resourceTableHtml(title, pairs, note) {
  return `
    <div class="section-label">${escapeHtml(title)}</div>
    <div class="table-scroll">
      <table class="data">
        <thead><tr><th>quantity</th><th>value</th><th>tag</th></tr></thead>
        <tbody>${quantityRows(pairs)}</tbody>
      </table>
    </div>
    ${note ? `<div class="method-note">${escapeHtml(note)}</div>` : ""}
  `;
}

/** Verdicts read "simulable in app" rather than "simulable_in_app". */
function readableVerdict(verdict) {
  const text = String(verdict || "").trim();
  if (!text) return "—";
  return text.replace(/_/g, " ");
}

function renderResources(body, status, el) {
  const tags = body.tags || {};
  const feas = body.feasibility || {};
  const subject = body.subject || {};
  const active = body.active_space;

  const pairs = [
    ["spatial orbitals", body.n_spatial_orbitals, tags.n_spatial_orbitals],
    ["electrons", body.n_electrons, tags.n_electrons],
    ["full qubits", body.n_qubits_full, tags.n_qubits_full],
    ["reduced qubits", body.n_qubits_after_parity_reduction, tags.n_qubits_after_parity_reduction],
    ["statevector memory", body.statevector_memory_human, tags.statevector_memory_human],
    ["verdict", readableVerdict(feas.verdict), tags["feasibility.verdict"]],
  ];

  /* The active-space block is opt-in and describes a smaller problem the caller
     cut down by hand, so every row in it is an estimate by definition. */
  const activePairs = active
    ? [
        ["active orbitals", active.active_orbitals, "estimated"],
        ["active qubits", active.n_qubits, "estimated"],
        ["reduced qubits", active.n_qubits_after_parity_reduction, "estimated"],
        ["statevector memory", active.statevector_memory_human, "estimated"],
        ["verdict", readableVerdict((active.feasibility || {}).verdict), "estimated"],
      ]
    : [];

  setStatus(status, "ok", "counted by formula · " + readableVerdict(feas.verdict));
  el.innerHTML = `
    <div class="computed-note">status: ${escapeHtml(body.status)}</div>
    ${resourceTableHtml(`${subject.formula || ""} · ${body.basis || ""}`, pairs, feas.reason)}
    <p class="muted">${escapeHtml(RESOURCE_FIXED_NOTE)}</p>
    ${active ? resourceTableHtml(`active space (k = ${active.active_orbitals} orbitals)`, activePairs, active.note) : ""}
    ${jsonBlock(body)}
  `;
}

async function runResources() {
  const formulaEl = $("#resource-formula");
  const activeEl = $("#resource-active-orbitals");
  const status = $("#resource-status");
  const detail = $("#resource-detail");
  const el = $("#resource-result");
  el.hidden = false;
  if (detail) detail.innerHTML = "";

  const formula = formulaEl && formulaEl.value ? formulaEl.value.trim() : "";
  const activeText = activeEl && activeEl.value ? String(activeEl.value).trim() : "";

  if (!formula) {
    setStatus(status, "err", "enter a formula");
    el.innerHTML = errorNote({
      code: "missing_formula",
      message: "Enter a formula such as H2O or C9H8O4 — the estimator counts from the formula, so it cannot guess one.",
    });
    return;
  }

  // The active-orbital count is optional and only sent when it is a positive
  // whole number, so an empty box means "full molecule" rather than k = 0.
  const params = new URLSearchParams({ formula });
  const active = parseInt(activeText, 10);
  if (activeText !== "" && isFinite(active) && active > 0) {
    params.set("active_orbitals", String(active));
  }

  setStatus(status, "placeholder", "counting…");
  try {
    const { body } = await apiFetch("/quantum/resources?" + params.toString());
    renderResources(body, status, el);
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

/* -------------------------------------------------------------- Molecule --- */

const molView = {
  viewer: null,
  spinning: false,
  structure3d: null,
  depiction2d: null,
  ensemble: null,
  conformers: [],
  activeConformer: 0,
  style: "ballandstick",
  labels: false,
  is3d: false,
  methodCode: null,
  methodDetail: null,
  elementLabels: [],
};

async function runSearch() {
  const q = $("#search-q").value;
  const status = $("#search-status");
  const el = $("#search-result");
  if (!q || !q.trim()) {
    setStatus(status, "err", "enter a query");
    return;
  }
  el.hidden = false;
  setStatus(status, "placeholder", "searching PubChem…");

  try {
    const { body } = await apiFetch(`/molecule/search?q=${encodeURIComponent(q)}`);
    setStatus(status, body.status === "computed" ? "ok" : "placeholder", body.status === "computed" ? `PubChem CID ${body.pubchem_cid}` : shortState(body));
    const smiles = body.smiles || "";
    el.innerHTML = (
      (body.status === "computed" ? `<div class="computed-note">status: ${body.status}</div>` : "") +
      (body.status === "placeholder" ? placeholderNote(body.note) : "") +
      `<dl>
        <dt>query</dt><dd>${escapeHtml(body.query)}</dd>
        <dt>CID</dt><dd>${body.pubchem_cid ?? "—"}</dd>
        <dt>SMILES</dt><dd>${escapeHtml(smiles || "—")}</dd>
        <dt>canonical</dt><dd>${escapeHtml(body.canonical_smiles || "—")}</dd>
        <dt>MW</dt><dd>${body.molecular_weight ?? "—"}</dd>
        <dt>formula</dt><dd>${escapeHtml(body.molecular_formula || "—")}</dd>
        <dt>IUPAC</dt><dd>${escapeHtml(body.iupac_name || "—")}</dd>
        <dt>3D SDF</dt><dd>${body.sdf_bytes && body.sdf_bytes.length ? `loaded (${body.sdf_dimension || "2D"}, ${body.sdf_bytes.length} bytes)` : "not loaded"}</dd>
      </dl>
      ${jsonBlock(body)}
      ${body.status === "computed" && smiles ? `<div class="row" style="margin-top:8px">
        <button id="use-smiles-btn" class="btn" type="button">Simulate this SMILES</button>
        ${body.sdf_bytes && body.sdf_bytes.length && body.sdf_dimension === "3D"
          ? `<button id="view-3d-btn" class="btn btn-primary" type="button">Load PubChem 3D conformer</button>`
          : (body.sdf_bytes && body.sdf_bytes.length
            ? `<span class="muted">PubChem returned a ${escapeHtml(body.sdf_dimension || "2D")} record — no 3D conformer to load.</span>`
            : "")}
      </div>` : ""}
    `);
    const useBtn = $("#use-smiles-btn");
    if (useBtn) useBtn.addEventListener("click", () => {
      $("#props-smiles").value = smiles;
      runMoleculeSimulation();
    });
    const btn = $("#view-3d-btn");
    if (btn) btn.addEventListener("click", () => loadMol3D(body));
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

async function runProperties() {
  const smiles = $("#props-smiles").value;
  const status = $("#props-status");
  const el = $("#props-result");
  if (!smiles || !smiles.trim()) {
    setStatus(status, "err", "enter a SMILES");
    return;
  }
  el.hidden = false;
  setStatus(status, "placeholder", "computing descriptors…");

  try {
    const { body } = await apiFetch("/molecule/properties", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles }),
    });

    setStatus(status, body.status === "computed" ? "ok" : "placeholder", body.status === "computed" ? "descriptors computed" : shortState(body));
    const lip = body.lipinski_rule_of_five;
    el.innerHTML = (
      (body.status === "computed" ? `<div class="computed-note">status: ${body.status}</div>` : "") +
      (body.status === "placeholder" ? placeholderNote(body.note) : "") +
      `<dl>
        <dt>SMILES</dt><dd>${escapeHtml(body.smiles)}</dd>
        <dt>MW</dt><dd>${body.molecular_weight ?? "—"}</dd>
        <dt>logP</dt><dd>${body.logp ?? "—"}</dd>
        <dt>HBD</dt><dd>${body.h_bond_donors ?? "—"}</dd>
        <dt>HBA</dt><dd>${body.h_bond_acceptors ?? "—"}</dd>
        <dt>rotatable</dt><dd>${body.rotatable_bonds ?? "—"}</dd>
        <dt>exact mass</dt><dd>${body.exact_mass ?? "—"}</dd>
        <dt>TPSA</dt><dd>${body.tpsa ?? "—"}</dd>
        <dt>formula</dt><dd>${escapeHtml(body.formula || "—")}</dd>
        <dt>Lipinski</dt><dd>${lip && typeof lip.pass === "boolean" ? (lip.pass ? "pass (" + lip.violations + " violations)" : "fail (" + lip.violations + " violations)") : "—"}</dd>
      </dl>
      ${jsonBlock(body)}
    `);
    // The formula is already on screen, so it costs one request to say what the
    // Quantum Lab would make of it. The line reports itself; a failure to count
    // (an element outside the estimator's table) is stated, not swallowed.
    if (body.status === "computed" && body.formula) {
      el.appendChild(quantumFeasibilityLine(body.formula));
    }
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

/**
 * One line for the Molecule Explorer: what the Quantum Lab's STO-3G resource
 * estimator says about a molecule whose formula is already shown. Nothing is
 * simulated for it and no number is invented here -- the counts come from the
 * same endpoint the Resource estimator panel uses.
 */
function quantumFeasibilityLine(formula) {
  const line = document.createElement("div");
  line.className = "method-note";
  line.id = "props-quantum-feasibility";
  line.textContent = "Quantum Lab feasibility: counting…";
  (async () => {
    try {
      const { body } = await apiFetch("/quantum/resources?formula=" + encodeURIComponent(formula));
      const feas = body.feasibility || {};
      line.textContent =
        `Quantum Lab feasibility: ${body.n_spatial_orbitals} spatial orbitals, ` +
        `${body.n_qubits_full} full / ${body.n_qubits_after_parity_reduction} reduced qubits — ` +
        `${readableVerdict(feas.verdict)} (STO-3G, counted by formula; no simulation was run).`;
    } catch (e) {
      const err = e && e.error ? e.error : e;
      line.textContent =
        `Quantum Lab feasibility: not available (${err.code || "error"}: ${err.message || ""}).`;
    }
  })();
  return line;
}

function verdictClass(verdict) {
  const v = String(verdict || "").toLowerCase();
  if (["pass", "likely", "high", "drug-like", "balanced"].includes(v)) return "good";
  if (["moderate", "unlikely", "very hydrophilic", "very lipophilic", "outside drug-like"].includes(v)) return "caution";
  if (["fail", "low", "adjacent"].includes(v)) return "bad";
  return "";
}

function findingsHtml(findings) {
  if (!Array.isArray(findings) || !findings.length) return "";
  return `
    <div class="section-label">Predicted behaviour (sourced rules)</div>
    <div class="findings">
      ${findings.map((f) => `
        <div class="finding ${verdictClass(f.verdict)}">
          <div class="finding-head">
            <strong>${escapeHtml(f.name)}</strong>
            <span class="finding-verdict">${escapeHtml(f.verdict || "")}</span>
          </div>
          ${f.value ? `<div class="finding-value">${escapeHtml(f.value)}</div>` : ""}
          ${f.rationale ? `<div class="finding-rationale">${escapeHtml(f.rationale)}</div>` : ""}
          ${f.source ? `<div class="finding-source">source: ${escapeHtml(f.source)}</div>` : ""}
        </div>
      `).join("")}
    </div>
  `;
}

function stepsHtml(steps, idPrefix) {
  if (!Array.isArray(steps) || !steps.length) return "";
  return `
    <div class="section-label">What the simulation did</div>
    <div class="steps" id="${idPrefix}-steps">
      ${steps.map((s, i) => `
        <div class="step" data-step-index="${i}" data-step="${escapeHtml(s.step || "")}">
          <div class="step-index">${i + 1}</div>
          <div class="step-body">
            <div class="step-title">
              ${escapeHtml(s.title || s.step || "")}
              ${s.status ? `<span class="chip ${escapeHtml(s.status)}">${escapeHtml(s.status)}</span>` : ""}
            </div>
            ${s.detail ? `<div class="step-detail">${escapeHtml(s.detail)}</div>` : ""}
            ${s.value ? `<div class="step-detail"><strong>${escapeHtml(s.value)}</strong></div>` : ""}
          </div>
        </div>
      `).join("")}
    </div>
  `;
}

function ensembleHtml(ensemble) {
  if (!ensemble || !ensemble.available) {
    return ensemble && ensemble.note ? `<div class="method-note">${escapeHtml(ensemble.note)}</div>` : "";
  }
  const rows = (ensemble.conformers || []).map((c) => `
    <tr data-conformer-index="${c.index}">
      <td class="numeric">${c.index}</td>
      <td class="numeric">${c.relative_energy_kcal_per_mol == null ? "—" : num(c.relative_energy_kcal_per_mol, 3)}</td>
      <td class="numeric">${c.population == null ? "—" : (100 * c.population).toFixed(1) + " %"}</td>
      <td>${c.minimised ? "yes" : "no"}</td>
    </tr>
  `).join("");
  return `
    <div class="section-label">Conformer ensemble (Boltzmann-weighted at ${num(ensemble.temperature_k, 2)} K)</div>
    <div class="table-scroll">
      <table class="data">
        <thead><tr><th>#</th><th>rel. energy (kcal/mol)</th><th>population (estimate)</th><th>minimised</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
    <div class="method-note">${escapeHtml(ensemble.method || "")}</div>
  `;
}

async function runMoleculeSimulation() {
  const smiles = $("#props-smiles").value;
  const status = $("#sim-status");
  const el = $("#sim-result");
  if (!smiles || !smiles.trim()) {
    setStatus(status, "err", "enter a SMILES");
    return;
  }
  const conformers = parseInt($("#sim-conformers").value, 10) || 8;
  el.hidden = false;
  setStatus(status, "placeholder", "embedding and minimising conformers…");

  try {
    const { body } = await apiFetch("/molecule/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles, num_conformers: conformers }),
    });

    if (body.status !== "computed") {
      setStatus(status, "placeholder", shortState(body));
      el.innerHTML = placeholderNote(body.note) + jsonBlock(body);
      return;
    }

    const p = body.properties || {};
    const st = body.structure || {};
    const ens = body.ensemble || {};
    const detail = [
      { k: "heavy atoms", v: p.heavy_atom_count },
      /* The 2D fallback ensemble carries no conformer_count at all, and
         String(undefined) would print the literal word "undefined". Say the
         count when there is one; otherwise name the situation instead. */
      { k: "conformers found", v: ens.conformer_count == null ? (st.is_3d ? "\u2014" : "not applicable (2D fallback)") : ens.conformer_count },
      { k: "energy span", v: ens.energy_span_kcal_per_mol == null ? "—" : num(ens.energy_span_kcal_per_mol, 2), sub: "kcal/mol (relative)" },
      { k: "MW", v: num(p.molecular_weight, 2) },
      { k: "logP", v: p.logp == null ? "—" : num(p.logp, 2) },
      { k: "TPSA", v: num(p.tpsa, 2) },
      { k: "HBD / HBA", v: `${p.h_bond_donors ?? "—"} / ${p.h_bond_acceptors ?? "—"}` },
      { k: "rotatable bonds", v: p.rotatable_bonds ?? "—" },
      { k: "rings", v: `${p.ring_count ?? "—"} (${p.aromatic_ring_count ?? "—"} aromatic)` },
    ];

    setStatus(status, "ok", st.is_3d ? "3D simulation complete" : "simulation complete (2D fallback)");
    el.innerHTML = `
      <div class="computed-note">status: ${escapeHtml(body.status)}</div>
      <div class="section-label">Result</div>
      <div class="metric-grid">
        ${detail.map((d) => `
          <div class="metric">
            <div class="k">${escapeHtml(d.k)}</div>
            <div class="v">${escapeHtml(d.v == null ? "\u2014" : String(d.v))}</div>
            ${d.sub ? `<div class="sub">${escapeHtml(d.sub)}</div>` : ""}
          </div>
        `).join("")}
      </div>
      <p class="method-note">${escapeHtml(body.interpretation_note || "")}</p>
      <div class="method-note"><strong>canonical SMILES:</strong> ${escapeHtml(body.canonical_smiles || body.smiles)}</div>
      ${findingsHtml(body.findings)}
      ${ensembleHtml(body.ensemble)}
      ${stepsHtml(body.simulation_steps, "sim")}
      ${jsonBlock(body)}
    `;

    molView.structure3d = st;
    molView.depiction2d = body.depiction || null;
    molView.ensemble = body.ensemble || null;
    showMolecule3D(st, body.method, body.method_detail);
    $("#mol-reset").disabled = false;
    $("#mol-2d").disabled = !(body.depiction && body.depiction.available);
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

/* ---------------------------------------------------------- Molecule 3D --- */

/** How far the loaded coordinates actually reach in z, in ångström. */
function zSpan(viewer) {
  let min_z = 0;
  let max_z = 0;
  let any = false;
  allAtoms(viewer).forEach((a) => {
    if (!isFinite(a.z)) return;
    any = true;
    if (a.z < min_z) min_z = a.z;
    if (a.z > max_z) max_z = a.z;
  });
  if (!any) return 0;
  return max_z - min_z;
}

/** Depth check on a molfile's own coordinates, so a flat record is never drawn as 3D. */
function molblockHasDepth(molblock, threshold = 0.25) {
  const lines = String(molblock || "").split(/\r?\n/);
  if (lines.length < 5) return false;
  const n = parseInt(String(lines[3]).slice(0, 3), 10);
  if (!isFinite(n) || n <= 0) return false;
  const zs = [];
  for (let i = 0; i < n && 4 + i < lines.length; i++) {
    const z = parseFloat(String(lines[4 + i]).slice(20, 30));
    if (isFinite(z)) zs.push(z);
  }
  if (!zs.length) return false;
  const min_z = Math.min(...zs);
  const max_z = Math.max(...zs);
  // Use z span (max z - min z), not max |z|, so a flat record whose z values
  // are numerically non-zero but identical is still reported as flat.
  return (max_z - min_z) > threshold;
}

/**
 * The badge answers "am I looking at a 3D structure or a drawing?" before
 * anything else. It follows the method code the API reports; when the API
 * reports none (a PubChem record, say) the depth is read off the coordinates
 * themselves, so a flat record is never presented as a 3D structure.
 */
function setMethodBadge(methodCode, methodDetail, viewer) {
  const badge = $("#mol-method-badge");
  if (!badge) return;
  let kind = "none";
  let text = "No structure is available in this environment";

  if (methodCode && methodCode.indexOf("3d_") === 0) {
    kind = "three-d";
    // Short label on the badge; the full sentence stays as the tooltip, because
    // a three-line badge covering the molecule is not a badge.
    text = methodCode === "3d_mmff94" ? "3D structure · MMFF94-minimised"
      : methodCode === "3d_uff" ? "3D structure · UFF-minimised"
      : methodCode === "3d_unminimised" ? "3D structure · not energy-minimised"
      : "3D structure";
  } else if (methodCode === "2d_fallback") {
    kind = "two-d";
    text = "Flat 2D layout, not a 3D structure";
  } else if (!methodCode && viewer) {
    const span = zSpan(viewer);
    if (span > 0.25) {
      kind = "three-d";
      text = "3D coordinates · " + num(span, 2) + " Å deep";
    } else {
      kind = "two-d";
      text = "Flat 2D layout, not a 3D structure";
    }
  }
  badge.textContent = text;
  badge.title = methodDetail || text;
  badge.dataset.kind = kind;
  badge.hidden = false;
}

function applyMoleculeStyle() {
  const viewer = molView.viewer;
  if (!viewer) return;
  try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
  if (!molView.is3d) {
    // A flat layout drawn as ball-and-stick would imply depth it does not have.
    viewer.setStyle({}, { stick: { radius: 0.16, colorscheme: ELEMENT_SCHEME } });
  } else {
    viewer.setStyle({}, styleSpec(molView.style));
    addHydrogenStyle(viewer);
    if (molView.style === "surface") translucentSurface(viewer);
  }
  setElementLabels(viewer, molView.labels ? allAtoms(viewer) : null, molView.elementLabels);
}

function redrawMolecule() {
  if (!molView.viewer) return;
  const styleEl = $("#mol-style");
  const labelsEl = $("#mol-labels");
  if (styleEl) molView.style = styleEl.value;
  if (labelsEl) molView.labels = labelsEl.checked;
  applyMoleculeStyle();
  try { molView.viewer.render(); } catch (e) { /* ignore */ }
}

function drawMolecule(molblock, format, is3d) {
  const viewer = molView.viewer;
  if (!viewer) return false;
  viewer.removeAllModels();
  viewer.removeAllLabels();
  try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
  molView.elementLabels.length = 0;
  viewer.addModel(molblock, format || "mol");
  molView.is3d = !!is3d;
  applyMoleculeStyle();
  frameAndTilt(viewer);
  viewer.render();
  return true;
}

function viewerLabel(viewer, text) {
  try {
    const atoms = centroid(viewer, {});
    if (!atoms) return;
    viewer.addLabel(text, {
      position: atoms,
      fontColor: "white",
      backgroundColor: "#1c2337",
      backgroundOpacity: 0.75,
      fontSize: 11,
    });
    viewer.render();
  } catch (e) {
    /* labels are cosmetic; never fail a render because of one */
  }
}

/**
 * Fill (or empty) the conformer selector. Each entry carries its own molblock,
 * so stepping through the ensemble shows real, separately minimised geometries
 * rather than an animation of one structure.
 */
function setConformers(ensemble) {
  const sel = $("#mol-conformer");
  if (!sel) return;
  const confs = ensemble && ensemble.available && Array.isArray(ensemble.conformers)
    ? ensemble.conformers.filter((c) => c && c.molblock)
    : [];
  molView.conformers = confs;
  molView.activeConformer = 0;
  sel.innerHTML = "";
  if (!confs.length) {
    const opt = document.createElement("option");
    opt.textContent = "not available";
    sel.appendChild(opt);
    sel.disabled = true;
    return;
  }
  sel.disabled = confs.length < 2;
  confs.forEach((c, i) => {
    const opt = document.createElement("option");
    opt.value = String(i);
    opt.textContent = `${i + 1} of ${confs.length}` + conformerSuffix(c);
    sel.appendChild(opt);
  });
  sel.value = "0";
}

function conformerSuffix(c) {
  const rel = c && c.relative_energy_kcal_per_mol;
  const pop = c && c.population;
  return (rel == null ? "" : ` · ${num(rel, 2)} kcal/mol`) +
    (pop == null ? "" : ` · ${(100 * pop).toFixed(1)} %`);
}

function showConformer(i) {
  const conf = molView.conformers[i];
  if (!conf || !molView.viewer || !conf.molblock) return;
  molView.activeConformer = i;
  if (!drawMolecule(conf.molblock, "mol", true)) return;
  viewerLabel(molView.viewer, `conformer ${i + 1}` + conformerSuffix(conf));
  setMethodBadge(molView.methodCode, molView.methodDetail, molView.viewer);
  setStatus($("#viewer-status"), "ok", `conformer ${i + 1} of ${molView.conformers.length}`);
}

function showMolecule3D(structure, methodCode, methodDetail) {
  const status = $("#viewer-status");
  const hint = $("#viewer-hint");
  if (!structure || !structure.molblock) {
    setStatus(status, "placeholder", "no structure produced");
    hint.textContent = "The simulation did not produce a structure to display.";
    setMethodBadge(null, null, null);
    setConformers(null);
    return;
  }
  const el = $("#mol-viewer");
  if (!molView.viewer) {
    const viewer = createViewer(el);
    if (!viewer) {
      setStatus(status, "err", "3Dmol library not loaded");
      hint.textContent = "3Dmol.js did not load; check the network tab for 3Dmol-min.js.";
      return;
    }
    molView.viewer = viewer;
    molView.spinning = false;
  }
  try {
    const is3d = !!structure.is_3d;
    molView.methodCode = is3d ? (methodCode || null) : (methodCode || "2d_fallback");
    molView.methodDetail = methodDetail || structure.method || null;
    drawMolecule(structure.molblock, "mol", is3d);
    setMethodBadge(molView.methodCode, molView.methodDetail, molView.viewer);
    setConformers(molView.ensemble);
    const controls = $("#mol-controls");
    if (controls) controls.hidden = false;
    setStatus(status, "ok", is3d ? "3D structure loaded" : "2D depiction loaded");
    hint.textContent = is3d
      ? "rotate with the mouse · the starting tilt shows the depth"
      : "this is a flat depiction, not a 3D structure";
  } catch (e) {
    setStatus(status, "err", "3Dmol viewer unavailable: " + (e && e.message ? e.message : e));
    hint.textContent = "3Dmol.js did not initialize; check the console for details.";
  }
}

function show2DDepiction() {
  const status = $("#viewer-status");
  const dep = molView.depiction2d;
  if (!dep || !dep.available || !dep.molblock) {
    setStatus(status, "placeholder", "no 2D depiction available");
    return;
  }
  if (molView.spinning) toggleSpin(molView, false);
  if (!drawMolecule(dep.molblock, "mol", false)) {
    setStatus(status, "err", "viewer unavailable");
    return;
  }
  molView.methodCode = "2d_fallback";
  setMethodBadge("2d_fallback", null, null);
  setStatus(status, "ok", "2D depiction loaded");
  $("#viewer-hint").textContent = "RDKit 2D layout · this is a flat depiction, not a 3D structure";
}

function toggleSpin(state, force) {
  if (!state || !state.viewer) return;
  state.spinning = typeof force === "boolean" ? force : !state.spinning;
  try {
    spinViewer(state.viewer, state.spinning);
    state.viewer.render();
  } catch (e) {
    state.spinning = false;
  }
}

function resetMolView() {
  if (!molView.viewer) return;
  try {
    frameAndTilt(molView.viewer);
    molView.viewer.render();
  } catch (e) { /* ignore */ }
}

function loadMol3D(body) {
  if (!body || !body.sdf_bytes || !body.sdf_bytes.length) {
    setStatus($("#viewer-status"), "placeholder", "no 3D data");
    return;
  }
  const el = $("#mol-viewer");
  if (!molView.viewer) {
    const viewer = createViewer(el);
    if (!viewer) {
      setStatus($("#viewer-status"), "err", "3Dmol library not loaded");
      $("#viewer-hint").textContent = "3Dmol.js did not load; check console and network for 3Dmol-min.js.";
      return;
    }
    molView.viewer = viewer;
    molView.spinning = false;
  }
  try {
    // The record is drawn as 3D only if its own coordinates have depth;
    // a flat PubChem record is drawn as the flat thing it is.
    const text = toText(body.sdf_bytes);
    const is3d = body.sdf_dimension === '3D' && molblockHasDepth(text);
    const title = is3d ? "PubChem 3D conformer (real)" : "PubChem 3D conformer (flat)";
    molView.methodCode = null;
    molView.methodDetail = null;
    molView.ensemble = null;

    if (!drawMolecule(text, "sdf", is3d)) {
      setStatus($("#viewer-status"), "err", "3Dmol failed to draw the loaded record");
      $("#viewer-hint").textContent = title;
      $("#mol-reset").disabled = false;
      return;
    }

    setMethodBadge(null, null, molView.viewer);
    setConformers(null);
    const controls = $("#mol-controls");
    if (controls) controls.hidden = false;
    setStatus($("#viewer-status"), "ok", is3d ? "3D record loaded" : "record loaded (flat)");
    $("#viewer-hint").textContent = is3d
      ? "PubChem 3D conformer loaded · rotate/zoom with mouse"
      : "PubChem returned a flat record — nothing here is 3D";
    $("#mol-reset").disabled = false;
  } catch (e) {
    setStatus($("#viewer-status"), "err", "3Dmol viewer unavailable: " + (e && e.message ? e.message : e));
    $("#viewer-hint").textContent = "3Dmol.js did not initialize; check console for details.";
  }
}

/* --------------------------------------------------------------- Protein --- */

async function fetchProtein() {
  const id = $("#protein-id").value.trim().toUpperCase();
  const status = $("#protein-status");
  const el = $("#protein-result");
  if (!id) {
    setStatus(status, "err", "enter a PDB id");
    return;
  }
  el.hidden = false;
  setStatus(status, "placeholder", "fetching RCSB PDB…");

  try {
    const { body } = await apiFetch(`/protein/${encodeURIComponent(id)}`);
    setStatus(status, body.status === "computed" ? "ok" : "placeholder", body.status === "computed" ? "entry loaded" : shortState(body));
    el.innerHTML = (
      (body.status === "computed" ? `<div class="computed-note">status: ${body.status}</div>` : "") +
      (body.status === "placeholder" ? placeholderNote(body.note) : "") +
      `<dl>
        <dt>PDB id</dt><dd>${escapeHtml(body.pdb_id)}</dd>
        <dt>title</dt><dd>${escapeHtml(body.title || "—")}</dd>
        <dt>resolution</dt><dd>${body.resolution ?? "—"} Å</dd>
        <dt>organism</dt><dd>${escapeHtml(body.organism || "—")}</dd>
        <dt>chains</dt><dd>${body.chains && body.chains.length ? body.chains.join(", ") : "—"}</dd>
        <dt>ligands</dt><dd>${body.ligands && body.ligands.length ? body.ligands.join(", ") : "—"}</dd>
        <dt>structure file</dt><dd>${escapeHtml(body.structure_file_reference || "—")}</dd>
      </dl>
      ${body.status === "computed" && body.structure_file_bytes && body.structure_file_bytes.length ? `<div class="computed-note" style="margin-top:8px">structure file loaded — click “View structure”</div>
       <button id="protein-view-btn" class="btn btn-primary" type="button">View structure</button>` : ""}
      ${jsonBlock(body)}
    `);
    const btn = $("#protein-view-btn");
    if (btn) btn.addEventListener("click", () => loadProtein3D(body));
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    el.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

const proteinView = {
  viewer: null,
  spinning: false,
  data: null,
  style: "cartoon",
  labels: false,
  elementLabels: [],
};

function redrawProtein() {
  if (!proteinView.viewer) return;
  const styleEl = $("#protein-style");
  const labelsEl = $("#protein-labels");
  if (styleEl) proteinView.style = styleEl.value;
  if (labelsEl) proteinView.labels = labelsEl.checked;
  applyProteinStyle(proteinView.viewer, {
    style: proteinView.style,
    labels: proteinView.labels,
    labelStore: proteinView.elementLabels,
  });
  try { proteinView.viewer.render(); } catch (e) { /* ignore */ }
}

function loadProtein3D(body) {
  if (!body || !body.structure_file_bytes || !body.structure_file_bytes.length) {
    setStatus($("#protein-viewer-status"), "placeholder", "no structure file");
    return;
  }
  setStatus($("#protein-viewer-status"), "placeholder", "building 3Dmol viewer…");
  const el = $("#protein-viewer");
  if (!proteinView.viewer) {
    const viewer = createViewer(el);
    if (!viewer) {
      setStatus($("#protein-viewer-status"), "err", "3Dmol library not loaded");
      $("#protein-viewer-hint").textContent = "3Dmol.js did not load; check console and network for 3Dmol-min.js.";
      return;
    }
    proteinView.viewer = viewer;
    proteinView.spinning = false;
  }
  try {
    const viewer = proteinView.viewer;
    const ref = String(body.structure_file_reference || "").toLowerCase();
    const ext = ref.endsWith(".pdb") || ref.endsWith(".ent") ? "pdb" : "cif";
    viewer.removeAllModels();
    viewer.removeAllLabels();
    try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
    proteinView.elementLabels.length = 0;
    viewer.addModel(toText(body.structure_file_bytes), ext);
    proteinView.data = body;
    redrawProtein();
    frameAndTilt(viewer);
    viewer.render();
    const controls = $("#protein-controls");
    if (controls) controls.hidden = false;
    const chains = chainsOf(viewer);
    setStatus($("#protein-viewer-status"), "ok", "structure loaded");
    $("#protein-viewer-hint").textContent =
      `cartoon coloured by chain (${chains.join(", ") || "—"}) · bound ligands as ball-and-stick · rotate/zoom with mouse`;
  } catch (e) {
    setStatus($("#protein-viewer-status"), "err", "3Dmol viewer unavailable: " + (e && e.message ? e.message : e));
    $("#protein-viewer-hint").textContent = "3Dmol.js did not initialize; check console for details.";
  }
}

/* --------------------------------------------------------------- Complex --- */

const complexView = {
  viewer: null,
  spinning: false,
  data: null,
  playToken: 0,
  style: "cartoon",
  labels: false,
  elementLabels: [],
};

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

const complexesState = { all: [], category: "all", query: "" };

function complexMatches(c) {
  const category = complexesState.category;
  if (category === "viral" && c.category !== "viral") return false;
  if (category === "human" && !String(c.category || "").startsWith("human")) return false;
  const query = complexesState.query.trim().toLowerCase();
  if (!query) return true;
  const haystack = `${c.display_name || ""} ${c.protein_pdb || ""} ${c.ligand_name || ""} ${c.ligand_code || ""}`.toLowerCase();
  return haystack.includes(query);
}

/** "Curated mechanism" vs "Contacts only" -- the tier badge the data declares. */
function tierBadgeHtml(tier) {
  return tier === "curated"
    ? '<span class="tier-badge curated">Curated mechanism</span>'
    : '<span class="tier-badge contacts_only">Contacts only</span>';
}

function renderComplexList() {
  const status = $("#complexes-status");
  const list = $("#complexes-list");
  if (!list) return;
  const total = complexesState.all.length;
  if (!total) {
    setStatus(status, "placeholder", "no curated complexes");
    list.innerHTML = `<div class="placeholder-note">No curated complexes available.</div>`;
    return;
  }
  const items = complexesState.all.filter(complexMatches);
  if (!items.length) {
    setStatus(status, "placeholder", `0 of ${total} complexes match`);
    list.innerHTML = '<div class="empty-filter-note">No curated complex matches this filter.</div>';
    return;
  }
  setStatus(status, "ok", `${items.length} of ${total} complexes`);
  list.innerHTML = items.map((c) => `
    <div class="list-item" data-complex-id="${escapeHtml(c.id || "")}">
      <div class="list-item-head">
        <h4>${escapeHtml(c.display_name || c.id)}</h4>
        ${tierBadgeHtml(c.tier)}
      </div>
      <p>${escapeHtml(c.protein_pdb || "")} · ${escapeHtml(c.ligand_name || "")}${c.ligand_code ? ` (${escapeHtml(c.ligand_code)})` : ""}${c.resolution_angstrom != null ? ` · ${num(c.resolution_angstrom, 2)} Å` : ""}</p>
      <p>${escapeHtml(c.target_class || "")}</p>
    </div>
  `).join("");
}

async function loadComplexes() {
  const status = $("#complexes-status");
  const list = $("#complexes-list");
  setStatus(status, "placeholder", "loading manifest…");
  try {
    const { body } = await apiFetch("/complexes/list");
    complexesState.all = Array.isArray(body.complexes) ? body.complexes : [];
    renderComplexList();
  } catch (e) {
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    list.innerHTML = errorNote(e && e.error ? e.error : e);
  }
}

function onComplexClick(e) {
  const item = e.target.closest(".list-item");
  if (!item) return;
  const id = item.dataset.complexId;
  $$(".list-item").forEach((i) => i.classList.remove("selected"));
  item.classList.add("selected");
  loadComplexDetail(id);
}

/** Per-type contact rules, each labelled as a measurement or an inference. */
function contactRulesHtml(rules, labels) {
  if (!rules || typeof rules !== "object") return "";
  const keys = Object.keys(rules);
  if (!keys.length) return "";
  const labelFor = (k) => (labels && labels[k] ? labels[k] : "measured geometry");
  return `<dl>${keys
    .map(
      (k) =>
        `<dt>${escapeHtml(String(k).replace(/_/g, " "))} · ${escapeHtml(labelFor(k))}</dt>` +
        `<dd>${escapeHtml(rules[k])}</dd>`,
    )
    .join("")}</dl>`;
}

function contactTableHtml(residues, idPrefix) {
  if (!Array.isArray(residues) || !residues.length) return "";
  const rows = residues.map((r, i) => `
    <tr data-contact-index="${i}" data-res-label="${escapeHtml(r.residue_label)}">
      <td>${escapeHtml(r.residue_label)}</td>
      <td class="numeric">${num(r.min_distance_angstrom, 2)}</td>
      <td class="numeric">${r.contact_count ?? "—"}</td>
      <td>${(r.types || []).map((t) => `<span class="type-tag ${escapeHtml(t)}">${escapeHtml(String(t).replace(/_/g, " "))}</span>`).join("")}</td>
    </tr>
  `).join("");
  return `
    <div class="section-label">Measured contacts (${residues.length} residues)</div>
    <div class="table-scroll">
      <table class="data" id="${idPrefix}-contact-table">
        <thead><tr><th>residue</th><th>min distance (Å)</th><th>contacts</th><th>interaction types</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  `;
}

function catalyticHtml(catalytic) {
  if (!Array.isArray(catalytic) || !catalytic.length) return "";
  return `
    <div class="section-label">Catalytic residues (literature) vs computed contact</div>
    <div class="catalytic-list">
      ${catalytic.map((c) => `
        <div class="catalytic ${c.in_contact ? "hit" : "miss"}">
          <span>
            <strong>${escapeHtml(c.label)}</strong>
            ${c.role ? `<span class="muted"> · ${escapeHtml(c.role)}</span>` : ""}
          </span>
          ${c.in_contact
            ? `<span class="dist">${num(c.min_distance_angstrom, 2)} Å in contact</span>`
            : `<span class="muted">not in contact</span>`}
        </div>
      `).join("")}
    </div>
  `;
}

function citationsHtml(citations) {
  if (!Array.isArray(citations) || !citations.length) return "";
  return `
    <div class="section-label">Sources (LITERATURE)</div>
    <ul class="citation-list">
      ${citations
        .map(
          (c) =>
            `<li>${escapeHtml(c.authors || "")} <em>${escapeHtml(c.title || "")}</em> ${escapeHtml(
              c.journal || "",
            )} ${escapeHtml(String(c.year || ""))}. ${
              c.doi
                ? `DOI <a href="https://doi.org/${encodeURIComponent(c.doi)}" target="_blank" rel="noopener">${escapeHtml(c.doi)}</a>`
                : ""
            }${c.pmid ? ` · PMID <a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(c.pmid)}/" target="_blank" rel="noopener">${escapeHtml(c.pmid)}</a>` : ""}</li>`,
        )
        .join("")}
    </ul>
  `;
}

function complexDetailHtml(body, id) {
  const mech = body.mechanism || {};
  // Only a curated entry gets a mechanism narrative; a contacts-only or explored
  // entry is measured geometry and nothing more.
  const curated = body.tier === "curated" && body.curation !== "none";
  const counts = body.contact_summary || {};
  const catalytic = Array.isArray(mech.catalytic_residues) ? mech.catalytic_residues : [];

  const metrics = [
    { k: "ligand", v: body.ligand_name || body.ligand_code, sub: `${body.ligand_atom_count ?? "—"} atoms` },
    { k: "binding residues", v: body.binding_site_residue_count ?? "—" },
    { k: "hydrogen bonds", v: counts.hydrogen_bond ?? "—", sub: "geometric heuristic · upper bound" },
    { k: "salt bridges", v: counts.salt_bridge ?? "—", sub: "geometric heuristic · upper bound" },
    { k: "hydrophobic", v: counts.hydrophobic ?? "—", sub: "measured geometry" },
    curated && catalytic.length
      ? {
          k: "catalytic blocked",
          v: `${mech.catalytic_residues_occluded ?? "—"} / ${mech.catalytic_residues_curated ?? "—"}`,
        }
      : { k: "catalytic blocked", v: "not curated", sub: "no literature-verified residue list" },
  ];
  if (body.resolution_angstrom != null) {
    metrics.push({ k: "resolution", v: `${num(body.resolution_angstrom, 2)} Å` });
  }

  // A donor/acceptor pair can have ideal geometry and still be excluded, purely
  // because the distance is over the cutoff. Say so, instead of showing a bare 0.
  const nearMiss = Array.isArray(body.near_miss_hydrogen_bonds) ? body.near_miss_hydrogen_bonds : [];
  const nearMissHtml = nearMiss.length
    ? `<p><strong>Not counted as hydrogen bonds</strong> — donor/acceptor geometry is right, but the distance is over 3.5 Å: ${nearMiss
        .map(
          (n) =>
            `${escapeHtml(n.protein_residue)} · ligand ${escapeHtml(n.ligand_atom)} · ${num(n.distance_angstrom, 2)} Å · ${n.angle_deg == null ? "angle n/a" : num(n.angle_deg, 1) + "\u00b0"}`,
        )
        .join("; ")}</p>`
    : "";

  // Every ligand N/O is treated as both a possible donor and a possible acceptor,
  // because the deposited file carries no hydrogens and no protonation states. The
  // set of real donors/acceptors is a subset of that, so the count can only over-state
  // the truth. Call it an upper bound rather than implying an exact count.
  const boundCount = counts.hydrogen_bond ?? 0;
  const boundCountSalt = counts.salt_bridge ?? 0;
  const upperBoundHtml =
    `<p><strong>Both counts above are upper bounds.</strong> ` +
    `Each ligand N/O is treated as a possible donor and a possible acceptor, because the ` +
    `deposited coordinates contain no hydrogens and no protonation states. Real donors and ` +
    `acceptors are a subset of that, so some of the ${escapeHtml(String(boundCount))} hydrogen-bond ` +
    `and ${escapeHtml(String(boundCountSalt))} salt-bridge residues may not form that interaction ` +
    `in reality. A salt bridge is checked before a hydrogen bond, so a residue can appear in both ` +
    `counts: the two lists overlap and must not be summed.</p>`;

  const rulesHtml = contactRulesHtml(body.contact_rules, body.contact_labels);
  const methodHtml = body.method ? `<p>${escapeHtml(body.method)}</p>` : "";

  // The API caps each residue's atom_pairs list; say so rather than showing the
  // capped list next to complete counts as if it were the whole set.
  const truncNote = body.atom_pairs_truncated
    ? `<p class="method-note">Per-residue atom pairs are capped at 12 (atom_pairs_truncated is set). The contact counts above are the complete tally.</p>`
    : "";

  // A distinct table id keeps the curated and explored panels from sharing one
  // DOM id when both are on screen.
  const tablePrefix = String(id || "").startsWith("explore:") ? "explore" : "complex";

  const uncuratedHtml = curated
    ? ""
    : `<p class="method-note"><strong>${escapeHtml(
        body.curation_note || "Not curated: measured contacts only",
      )}</strong> There is no mechanism narrative and no catalytic-residue panel for this entry; everything below is measured geometry from the deposited coordinates.</p>`;

  return `
      <div class="computed-note">
        <span class="tier-badge ${curated ? "curated" : "uncurated"}">${
          curated ? "Curated mechanism" : "Not curated: measured contacts only"
        }</span>
        ${escapeHtml(body.display_name || id)}
      </div>
      ${uncuratedHtml}
      <div class="metric-grid" style="margin-top:10px">
        ${metrics
          .map(
            (m) => `
          <div class="metric">
            <div class="k">${escapeHtml(m.k)}</div>
            <div class="v">${escapeHtml(String(m.v))}</div>
            ${m.sub ? `<div class="sub">${escapeHtml(m.sub)}</div>` : ""}
          </div>
        `,
          )
          .join("")}
      </div>
      <p class="method-note">Every distance below is measured from the deposited crystal
      coordinates, in ångström, between heavy atoms only.</p>
      <details class="contact-rules">
        <summary>How these contacts are measured</summary>
        <div>
          ${methodHtml}
          ${rulesHtml}
          ${upperBoundHtml}
          ${nearMissHtml}
        </div>
      </details>
      ${curated ? citationsHtml(body.citations) : ""}
      ${curated ? catalyticHtml(catalytic) : ""}
      ${curated && body.catalytic_residues_note ? `<p class="method-note">${escapeHtml(body.catalytic_residues_note)}</p>` : ""}
      ${curated ? stepsHtml(mech.steps, "mech") : ""}
      ${contactTableHtml(body.interacting_residues, tablePrefix)}
      ${truncNote}
      ${jsonBlock(body)}
    `;
}

function renderComplexInto(el, body, id) {
  el.hidden = false;
  el.innerHTML = complexDetailHtml(body, id);
  complexView.data = body;
  setComplexControls(true);
  showComplex3D(body);
  const play = $("#complex-play");
  if (play) {
    const steps = body.mechanism && Array.isArray(body.mechanism.steps) ? body.mechanism.steps : [];
    play.disabled = !steps.length;
  }
}

async function loadComplexDetail(id) {
  const status = $("#complexes-status");
  const el = $("#complex-detail");
  el.hidden = false;
  if (!id) {
    el.innerHTML = '<p class="muted">Select a complex from the list.</p>';
    setStatus(status, "placeholder", "");
    return;
  }
  setStatus(status, "placeholder", `measuring contacts for ${id}…`);
  el.innerHTML = `<div class="method-note">Downloading the deposited coordinates and measuring contacts…</div>`;

  try {
    const { body } = await apiFetch(`/complexes/${encodeURIComponent(id)}/structure`);
    if (body.status !== "computed") {
      setStatus(status, "placeholder", shortState(body));
      el.innerHTML = placeholderNote(body.note) + jsonBlock(body);
      complexView.data = null;
      setComplexControls(false);
      return;
    }

    setStatus(status, "ok", `contacts measured for ${body.protein_pdb}`);
    renderComplexInto(el, body, id);
  } catch (e) {
    el.innerHTML = errorNote(e && e.error ? e.error : e);
    setStatus(status, "err", shortErrorState(e && e.error ? e.error : e));
    complexView.data = null;
    setComplexControls(false);
  }
}

function onExploreCandidateClick(e) {
  const button = e.target.closest("[data-explore-ligand]");
  if (!button) return;
  runExplore($("#explore-pdb").value, button.dataset.exploreLigand);
}

async function runExplore(pdbId, ligand) {
  const status = $("#explore-status");
  const candidates = $("#explore-candidates");
  const detail = $("#explore-detail");
  const id = String(pdbId || "").trim().toUpperCase();
  if (!/^[0-9][A-Za-z0-9]{3}$/.test(id)) {
    setStatus(status, "err", "enter a 4-character PDB ID");
    candidates.innerHTML = "";
    detail.hidden = true;
    return;
  }

  setStatus(status, "placeholder", `fetching ${id}…`);
  candidates.innerHTML = "";
  detail.hidden = true;
  try {
    const suffix = ligand ? `?ligand=${encodeURIComponent(ligand)}` : "";
    const { body } = await apiFetch(`/complexes/explore/${encodeURIComponent(id)}${suffix}`);

    if (body.status === "ligand_selection_required") {
      setStatus(status, "placeholder", `${body.ligand_candidates.length} drug-like ligands`);
      candidates.innerHTML = `
        <p class="method-note">${escapeHtml(body.note || "Choose a ligand.")}</p>
        <p class="method-note">${escapeHtml(body.ligand_rule || "")}</p>
        <div class="candidate-row">
          ${body.ligand_candidates
            .map(
              (c) =>
                `<button class="btn candidate-btn" type="button" data-explore-ligand="${escapeHtml(
                  c.ligand_code,
                )}">${escapeHtml(c.ligand_code)} · ${c.heavy_atom_count} heavy atoms${
                  c.chain ? ` · chain ${escapeHtml(c.chain)}` : ""
                }</button>`,
            )
            .join("")}
        </div>
      `;
      return;
    }

    if (body.status !== "computed") {
      setStatus(status, "placeholder", shortState(body));
      candidates.innerHTML = placeholderNote(body.note || "No contacts could be measured.");
      return;
    }

    setStatus(status, "ok", `${body.protein_pdb} · ${body.ligand_code} · contacts measured`);
    renderComplexInto(detail, body, body.complex_id);
  } catch (e) {
    const err = e && e.error ? e.error : e;
    setStatus(status, "err", shortErrorState(err));
    candidates.innerHTML = errorNote(err);
    detail.hidden = true;
  }
}

function setComplexControls(enabled) {
  $("#complex-reset").disabled = !enabled;
  $("#complex-play").disabled = !enabled;
  const controls = $("#complex-controls");
  if (controls) controls.hidden = !enabled;
}

function redrawComplex() {
  if (!complexView.viewer) return;
  const styleEl = $("#complex-style");
  const labelsEl = $("#complex-labels");
  if (styleEl) complexView.style = styleEl.value;
  if (labelsEl) complexView.labels = labelsEl.checked;
  complexStyles(complexView.viewer);
  try { complexView.viewer.render(); } catch (e) { /* ignore */ }
}

/**
 * Same style vocabulary as the other two viewers: the protein is a cartoon
 * coloured by chain, and the bound drug is always ball-and-stick on top of
 * whatever the protein is doing.
 */
function complexStyles(viewer) {
  try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
  applyProteinStyle(viewer, {
    style: complexView.style,
    labels: complexView.labels,
    labelStore: complexView.elementLabels,
  });
}

function showComplex3D(data) {
  const status = $("#complex-viewer-status");
  const hint = $("#complex-viewer-hint");
  const text = toText(data.structure_file_bytes);
  if (!text) {
    setStatus(status, "placeholder", "no structure file produced");
    return;
  }
  const el = $("#complex-viewer");
  if (!complexView.viewer) {
    const viewer = createViewer(el);
    if (!viewer) {
      setStatus(status, "err", "3Dmol library not loaded");
      hint.textContent = "3Dmol.js did not load; check the network tab for 3Dmol-min.js.";
      return;
    }
    complexView.viewer = viewer;
    complexView.spinning = false;
  }
  const viewer = complexView.viewer;
  try {
    viewer.removeAllModels();
    viewer.removeAllLabels();
    try { viewer.removeAllSurfaces(); } catch (e) { /* ignore */ }
    complexView.elementLabels.length = 0;
    viewer.addModel(text, "pdb");
    complexStyles(viewer);

    const ligCentre = centroid(viewer, { hetflag: true });
    if (ligCentre) {
      viewer.addLabel(`${data.ligand_name || "ligand"} (${data.ligand_code || "?"})`, {
        position: ligCentre,
        fontColor: "white",
        backgroundColor: "#1c2337",
        backgroundOpacity: 0.8,
        fontSize: 11,
      });
    }
    viewer.zoomTo({ hetflag: true });
    viewer.zoom(0.55);
    try {
      if (typeof viewer.rotate === "function") {
        viewer.rotate((TILT_Y_DEG * Math.PI) / 180, "y");
        viewer.rotate((TILT_X_DEG * Math.PI) / 180, "x");
      }
    } catch (e) { /* ignore */ }
    viewer.render();
    setStatus(status, "ok", `${data.protein_pdb} · protein cartoon + bound drug`);
    hint.textContent = "Protein as cartoon coloured by chain, drug as ball-and-stick · press “Play mechanism” to walk the measured contacts";
  } catch (e) {
    setStatus(status, "err", "3Dmol viewer unavailable: " + (e && e.message ? e.message : e));
    hint.textContent = "3Dmol.js did not initialize; check the console for details.";
  }
}

/** Parse "ASP A25" / "HIS A383" into a 3Dmol residue selection. */
function specFromResidueLabel(label) {
  const parts = String(label || "").trim().split(/\s+/);
  if (parts.length < 2) return null;
  const chainAndNumber = parts[1];
  const chain = chainAndNumber.replace(/[^A-Za-z]/g, "");
  const resi = chainAndNumber.replace(/[^0-9]/g, "");
  if (!resi) return null;
  const spec = { resi: resi };
  if (chain) spec.chain = chain;
  return spec;
}

async function playMechanism() {
  const data = complexView.data;
  const viewer = complexView.viewer;
  if (!data || !viewer) return;

  const token = ++complexView.playToken;
  const playBtn = $("#complex-play");
  const status = $("#complex-viewer-status");
  const residues = (data.interacting_residues || []).slice(0, 12);
  const catalytic = (data.mechanism && data.mechanism.catalytic_residues) || [];

  playBtn.disabled = true;
  playBtn.textContent = "Playing…";

  const markStep = (name) => {
    $$("#mech-steps .step").forEach((s) => {
      const active = s.dataset.step === name;
      s.classList.toggle("playing", active);
      if (active) s.scrollIntoView({ block: "nearest", behavior: "smooth" });
    });
  };
  const flashRow = (label) => {
    $$("#complex-contact-table tbody tr").forEach((tr) => {
      tr.classList.toggle("flash", tr.dataset.resLabel === label);
    });
  };
  const clearHighlights = () => {
    complexStyles(viewer);
    viewer.removeAllLabels();
  };

  try {
    spinViewer(viewer, true);

    // 1. The drug sits in the pocket.
    markStep("occupy");
    setStatus(status, "placeholder", "step 1 · drug in the binding pocket");
    const ligCentre = centroid(viewer, { hetflag: true });
    if (ligCentre) {
      viewer.addLabel(`${data.ligand_name || "drug"} occupying the pocket`, {
        position: ligCentre, fontColor: "white", backgroundColor: "#1c2337",
        backgroundOpacity: 0.8, fontSize: 11,
      });
    }
    viewer.render();
    await wait(1600);
    if (token !== complexView.playToken) return;

    // 2. Walk each measured contact, nearest first.
    for (let i = 0; i < residues.length; i++) {
      const r = residues[i];
      clearHighlights();
      const spec = residueSpec(r);
      viewer.addStyle(spec, HIGHLIGHT);
      const c = centroid(viewer, spec);
      if (c) {
        viewer.addLabel(`${r.residue_label} — ${num(r.min_distance_angstrom, 2)} Å`, {
          position: c, fontColor: "white", backgroundColor: "#1c2337",
          backgroundOpacity: 0.85, fontSize: 11,
        });
      }
      viewer.render();
      markStep("anchor");
      flashRow(r.residue_label);
      setStatus(status, "placeholder", `contact ${i + 1}/${residues.length} · ${r.residue_label} at ${num(r.min_distance_angstrom, 2)} Å`);
      await wait(900);
      if (token !== complexView.playToken) return;
    }

    // 3. Highlight the catalytic machinery the drug is sitting on.
    clearHighlights();
    viewer.spin(false);
    const hits = catalytic.filter((c) => c.in_contact);
    hits.forEach((c) => {
      const spec = specFromResidueLabel(c.match);
      if (!spec) return;
      viewer.addStyle(spec, CATALYTIC_HIGHLIGHT);
      const centre = centroid(viewer, spec);
      if (centre) {
        viewer.addLabel(`${c.label} — catalytic, ${num(c.min_distance_angstrom, 2)} Å`, {
          position: centre, fontColor: "white", backgroundColor: "#7a1f1f",
          backgroundOpacity: 0.85, fontSize: 11,
        });
      }
    });
    viewer.render();
    markStep("blockade");
    flashRow("");
    if (hits.length) {
      setStatus(status, "ok", `catalytic residues blocked: ${hits.map((c) => c.label).join(", ")}`);
    } else {
      setStatus(status, "placeholder", "no curated catalytic residue is inside the contact cutoff");
    }
    await wait(2200);
    if (token !== complexView.playToken) return;

    // 4. Consequence.
    clearHighlights();
    markStep("consequence");
    const centre = centroid(viewer, { hetflag: true });
    if (centre) {
      viewer.addLabel("active site blocked — substrate cannot be processed", {
        position: centre, fontColor: "white", backgroundColor: "#1c2337",
        backgroundOpacity: 0.85, fontSize: 11,
      });
    }
    viewer.render();
    setStatus(status, "ok", "mechanism complete");
  } finally {
    if (token === complexView.playToken) {
      playBtn.disabled = false;
      playBtn.textContent = "Play mechanism";
      try { viewer.spin(false); } catch (e) { /* ignore */ }
    }
  }
}

function resetComplexView() {
  const viewer = complexView.viewer;
  if (!viewer) return;
  complexView.playToken++;
  $("#complex-play").disabled = !complexView.data;
  $("#complex-play").textContent = "Play mechanism";
  try {
    viewer.spin(false);
    complexView.spinning = false;
    const rotateBox = $("#complex-rotate");
    if (rotateBox) rotateBox.checked = false;
    viewer.removeAllLabels();
    complexStyles(viewer);
    viewer.zoomTo({ hetflag: true });
    viewer.zoom(0.55);
    try {
      if (typeof viewer.rotate === "function") {
        viewer.rotate((TILT_Y_DEG * Math.PI) / 180, "y");
        viewer.rotate((TILT_X_DEG * Math.PI) / 180, "x");
      }
    } catch (e) { /* ignore */ }
    viewer.render();
    $$("#complex-contact-table tbody tr").forEach((tr) => tr.classList.remove("flash"));
    $$("#mech-steps .step").forEach((s) => s.classList.remove("playing"));
    setStatus($("#complex-viewer-status"), "ok", "view reset");
  } catch (e) { /* ignore */ }
}
