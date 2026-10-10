# tpt-optics task runner. Requires `just` (https://github.com/casey/just).
# All recipes delegate to `xtask` so CI and local dev stay identical.

default:
    @just --list

build:
    cargo run -p xtask -- build

test:
    cargo run -p xtask -- test

lint:
    cargo run -p xtask -- lint

fmt:
    cargo run -p xtask -- fmt

fmt-check:
    cargo run -p xtask -- fmt --check

deny:
    cargo run -p xtask -- deny

# Full local gate: format check, clippy, tests, license check, no_std build.
check:
    cargo run -p xtask -- check
