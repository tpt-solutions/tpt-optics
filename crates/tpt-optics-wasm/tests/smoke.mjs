// Smoke test for the Node build of tpt-optics-wasm.
//
// Build first:   wasm-pack build crates/tpt-optics-wasm --target nodejs --out-dir pkg-node
// Then run:      node crates/tpt-optics-wasm/tests/smoke.mjs

import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import assert from "node:assert/strict";

const require = createRequire(import.meta.url);
const here = dirname(fileURLToPath(import.meta.url));
const wasm = require(join(here, "..", "pkg-node", "tpt_optics_wasm.js"));

// N-BK7 at the He-Ne line: n = 1.5151 (Sellmeier model).
const n = wasm.glass_index("N-BK7", 632.8);
assert.ok(Math.abs(n - 1.5151) < 1e-3, `n = ${n}`);

// Unknown glasses are reported, not silently accepted.
assert.throws(() => wasm.glass_index("NOT-A-GLASS", 587.6));

// Singlet report: EFL close to the thick-lens value (about 97.6 mm).
const report = JSON.parse(wasm.singlet_report(100, -100, 5, 12.5, "N-BK7"));
assert.ok(report.efl > 0.09 && report.efl < 0.1, `EFL ${report.efl}`);
// The image plane is at the focus, so the back focal length from it is zero.
assert.ok(Math.abs(report.bfl) < 1e-6, `BFL ${report.bfl}`);
assert.equal(report.fields.length, 1);
// Image plane at focus: the uncorrected f/4 singlet's spot is a fraction of a
// millimetre. Defocused by the back focal length it is several millimetres.
assert.ok(report.fields[0].rms_spot < 5e-4, `spot ${report.fields[0].rms_spot}`);

console.log(
  `ok: n(N-BK7, 632.8 nm) = ${n.toFixed(5)}, EFL = ${(report.efl * 1e3).toFixed(2)} mm, ` +
    `on-axis RMS spot = ${(report.fields[0].rms_spot * 1e6).toFixed(2)} µm`,
);
