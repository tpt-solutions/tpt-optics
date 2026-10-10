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
