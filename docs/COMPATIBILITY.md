# File-format compatibility

What the `tpt-optics-io` readers and writers accept, what they reject, and what
they drop. Anything not listed as supported is an error, or is listed here as
ignored. A silent misreading is treated as a bug.

## Zemax `.zmx` (sequential, text)

| Keyword | Read | Written | Notes |
|---|---|---|---|
| `MODE SEQ` | yes | yes | Other modes (`MODE NSC`) are errors. |
| `NAME` (system or surface) | yes | yes | |
| `SURF` | yes | yes | One block per surface. |
| `CURV` | yes | yes | Radius in mm is converted to 1/m. The asphere terms after `CURV` are read in lens units. |
| `CONI` | yes | yes | |
| `DISZ` | yes | yes | `INFINITY` is supported. The object-space `DISZ` is kept on surface 0 but not used. |
| `GLAS` | yes | yes | A catalogue name, `AIR`, a plain index, or `MIRROR`. Unknown names are errors with a suggestion. |
| `DIAM` | yes | yes | Semi-diameter in mm. |
| `STOP` | yes | yes | Kept in the document. It does not change the trace. |
| `TYPE` | `STANDARD` only | no | Other types (for example `EVENASPH`) are errors. |
| `PARM` | zero values only | no | A non-zero value is an error, because it belongs to a non-standard surface type. |
| `FNUM`, `WAVL`, and other headers | ignored | no | |

Written numbers: curvature and thickness use the file's units. Asphere terms
use the shortest exact form, so they round-trip. Index hints use the shortest
exact form too.

A UTF-8 byte-order mark is accepted. UTF-16 input must be decoded first with
`tpt_optics::io::text::decode`, which accepts UTF-16 with a BOM.

Binary `.zmx` and non-sequential (`NSC`) files are out of scope.

## CODE V `.seq` (sequence)

| Command | Read | Written | Notes |
|---|---|---|---|
| `LEN NEW` | yes | yes | |
| `SO` | yes | yes | The object surface. |
| `S` | yes | yes | `radius thickness [glass]`. `INF` radius is flat, `INF` thickness is infinite. |
| `STO` | yes | yes | Sets the stop flag on the next surface. |
| `TITLE`, `RMS`, `GO`, `EFI`, `EFL`, `REF` | ignored | `GO` only | |
| Anything else | error | no | |

Radii and thicknesses are in millimetres in the file. `nan` and `inf` values
are errors, except the `INF` keyword.

## Glass catalogues (AGF)

| Record | Read | Notes |
|---|---|---|
| `NM` | yes | Name, and the dispersion formula number (the second field). |
| `CD` | yes | Coefficients may wrap over several lines. More than six is an error. |
| `ED` | yes | The dn/dT field is taken from position 5. Check that layout against a vendor file. |
| `TD` | ignored | Thermal data, not used. The temperature coefficient comes from `ED`. |
| `CC` | ignored | Concatenation, not used. |
| Formula 1 (Schott) | yes | `n² = A0 + A1 λ² + A2 λ⁻² + …` |
| Formula 2 (Sellmeier) | yes | `K1 L1 K2 L2 K3 L3` |
| Other formulas | error | Not read as Sellmeier. |

Vendor catalogue files are not bundled. The built-in reference catalogue has
a small set of glasses for tests and examples.

## Design files (TOML)

The CLI reads TOML design files. The schema is printed by `tpt-optics schema`
as JSON Schema, and the sample is `crates/tpt-optics-cli/examples/singlet.toml`.
Unknown operand kinds and unknown glasses are errors.

## Results

- JSON from `AnalysisReport::to_json` and the CLI `--json` flag. Non-finite
  numbers are written as `null`.
- CSV from `AnalysisReport::to_csv`, one row per field.
- SVG spot diagrams from `AnalysisReport::spot_svg`.
