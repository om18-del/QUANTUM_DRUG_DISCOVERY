#!/usr/bin/env python3
"""Fix the broken loadMol3D function in frontend/static/app.js.

The function got duplicated/corrupted during earlier edits. This script
replaces the entire broken section with a clean version.
"""
from pathlib import Path

APP_JS = Path(__file__).resolve().parents[1] / "frontend" / "static" / "app.js"

content = APP_JS.read_text(encoding="utf-8")

# Find and replace the corrupted loadMol3D function
old_func_start = "function loadMol3D(body) {"
old_func_end = "}\n\n/* ----"

start_idx = content.find(old_func_start)
if start_idx == -1:
    print("ERROR: Could not find loadMol3D function start")
    exit(1)

end_idx = content.find(old_func_end, start_idx)
if end_idx == -1:
    print("ERROR: Could not find loadMol3D function end")
    exit(1)

end_idx += len(old_func_end)

new_function = """function loadMol3D(body) {
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
    molView.methodCode = null;
    molView.methodDetail = null;
    molView.ensemble = null;

    if (!drawMolecule(text, "sdf", is3d)) {
      throw new Error('3Dmol failed to draw the loaded record');
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

/* ----"""

new_content = content[:start_idx] + new_function + content[end_idx:]
APP_JS.write_text(new_content, encoding="utf-8")
print("Fixed loadMol3D function in app.js")