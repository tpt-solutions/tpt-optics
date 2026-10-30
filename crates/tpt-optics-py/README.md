# tpt-optics-py

Python bindings for tpt-optics, built with PyO3. The module is `tpt_optics`.

This crate is not a workspace member. Its PyO3 build dependency,
`target-lexicon`, is Apache-2.0 with the LLVM exception. The workspace licence
policy does not allow that licence, so the crate builds on its own.

## Build

```sh
cd crates/tpt-optics-py
cargo build --release --features extension-module
```

On Windows, copy `target/release/tpt_optics_py.dll` to a directory on
`sys.path` as `tpt_optics.pyd`. On Linux and macOS, copy
`libtpt_optics_py.so` or `libtpt_optics_py.dylib` as `tpt_optics.so`. Or build
a wheel with `maturin build --release --features extension-module`.

## Use

```python
import json
import tpt_optics as tpt

tpt.glass_index("N-BK7", 632.8)             # refractive index
report = json.loads(tpt.singlet_report(100.0, -100.0, 5.0, 12.5, "N-BK7"))
report["efl"]                                # effective focal length, m
tpt.rank_achromat_pairs(100.0, ["bk7"], ["f2", "n-sf11"])
```

Lengths are in millimetres at this boundary. `singlet_report` returns the
analysis as JSON. An unknown glass raises `KeyError`, and an invalid
prescription raises `ValueError`.

## Test

With the extension importable, run `pytest tests/`.
