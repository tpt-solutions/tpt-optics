//! Repository automation for the tpt-optics workspace.
//!
//! Usage: `cargo run -p xtask -- <command>`
//!
//! Commands:
//! - `build` — build every workspace member
//! - `test`  — run the full test suite
//! - `fmt`   — format (`--check` to verify only)
//! - `lint`  — clippy with `-D warnings`
//! - `deny`  — `cargo deny check` (licenses, bans, sources)
//! - `nostd` — build the `no_std`-compatible crates without default features
//! - `check` — fmt --check + lint + test + deny + nostd (the CI gate)

use std::process::Command;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cmd = args.first().map_or("check", String::as_str);
    let rest = &args[1.min(args.len())..];
    let status = match cmd {
        "build" => cargo(&["build", "--workspace", "--all-targets"]),
        "test" => cargo(&["test", "--workspace"]),
        "fmt" => {
            let mut a = vec!["fmt", "--all"];
            if rest.iter().any(|r| r == "--check") {
                a.push("--check");
            }
            cargo(&a)
        }
        "lint" => cargo(&[
            "clippy",
            "--workspace",
            "--all-targets",
            "--",
            "-D",
            "warnings",
        ]),
        "deny" => cargo(&["deny", "check", "licenses", "sources"]),
        "nostd" => nostd(),
        "validation" => validation(),
        "check" => {
            type Step = (&'static str, fn() -> i32);
            let steps: [Step; 5] = [
                ("fmt --check", || cargo(&["fmt", "--all", "--check"])),
                ("lint", lint),
                ("test", || cargo(&["test", "--workspace"])),
                ("deny", || cargo(&["deny", "check", "licenses", "sources"])),
                ("nostd", nostd),
            ];
            let mut code = 0;
            for (name, step) in steps {
                println!("==> xtask check: {name}");
                if step() != 0 {
                    code = 1;
                    println!("==> step `{name}` FAILED");
                }
            }
            code
        }
        other => {
            eprintln!("unknown xtask command: {other}");
            2
        }
    };
    std::process::exit(status);
}

fn cargo(args: &[&str]) -> i32 {
    println!("==> cargo {}", args.join(" "));
    let status = Command::new("cargo").args(args).status();
    match status {
        Ok(s) if s.success() => 0,
        Ok(_) => 1,
        Err(e) => {
            eprintln!("failed to spawn cargo: {e}");
            127
        }
    }
}

fn lint() -> i32 {
    cargo(&[
        "clippy",
        "--workspace",
        "--all-targets",
        "--",
        "-D",
        "warnings",
    ])
}

/// Build the core crates with `--no-default-features` to prove they stay
/// `no_std`-compatible (any accidental `std` dependency fails the build).
fn nostd() -> i32 {
    let mut code = 0;
    for krate in [
        "tpt-optics-core",
        "tpt-optics-geometric",
        "tpt-optics-polarization",
    ] {
        if cargo(&["build", "-p", krate, "--no-default-features"]) != 0 {
            code = 1;
        }
    }
    code
}

/// Check that every test name cited in `docs/VALIDATION.md` exists in the
/// source. Names are the backticked snake_case identifiers in the page.
fn validation() -> i32 {
    let page = match std::fs::read_to_string("docs/VALIDATION.md") {
        Ok(p) => p,
        Err(e) => {
            eprintln!("cannot read docs/VALIDATION.md: {e}");
            return 1;
        }
    };
    let mut sources = String::new();
    collect_rs(std::path::Path::new("crates"), &mut sources);
    let mut missing = Vec::new();
    let mut checked = 0;
    for token in page.split('`').skip(1).step_by(2) {
        let is_name = token.contains('_')
            && token
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_');
        if !is_name {
            continue;
        }
        checked += 1;
        if !sources.contains(&format!("fn {token}(")) {
            missing.push(token.to_string());
        }
    }
    if missing.is_empty() {
        println!("validation: {checked} cited tests found");
        0
    } else {
        for name in &missing {
            eprintln!("validation: test '{name}' is cited in docs/VALIDATION.md but not found");
        }
        1
    }
}

fn collect_rs(dir: &std::path::Path, out: &mut String) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_rs(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            if let Ok(text) = std::fs::read_to_string(&path) {
                out.push_str(&text);
            }
        }
    }
}
