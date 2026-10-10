# Security Policy

## Supported versions

| Version | Supported |
|---------|-----------|
| 0.1.x   | ✅        |

## Reporting a vulnerability

Do **not** open a public issue for security problems.

Email **security@tpt.solutions** with:

- a description of the issue,
- steps or a minimal reproducer,
- affected versions and any known workarounds.

You will receive an acknowledgment within 5 business days. Please allow up to
90 days for a fix before public disclosure; we will credit reporters in the
release notes unless asked otherwise.

## Scope notes

tpt-optics is a numerical simulation library. It parses untrusted input in
exactly one place: the `tpt-optics-io` file-format parsers (Zemax `.zmx`,
CODE V `.seq`, Zemax AGF catalogs). Memory-safety bugs in those parsers are
in scope; because the workspace forbids `unsafe` code everywhere
(`unsafe_code = "forbid"`), memory safety violations would be compiler bugs —
panics, hangs, or unbounded memory growth triggered by malformed files are
the realistic class of issues. Denial-of-service via pathological but valid
inputs is treated as a regular bug, not a security vulnerability.
