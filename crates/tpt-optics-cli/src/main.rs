//! `tpt-optics` command-line tool.
//!
//! ```text
//! tpt-optics info
//! tpt-optics glass <NAME> [--nm <WAVELENGTH>]
//! tpt-optics analyze <DESIGN.toml> [--json | --csv | --svg <FIELD_INDEX>]
//! tpt-optics optimize <DESIGN.toml> [--write <OUT.toml>]
//! tpt-optics trace <DESIGN.toml> --height-mm <H> [--field-deg <F>]
//! tpt-optics run <DESIGN.toml> --out <DIR>
//! tpt-optics convert <IN> <OUT>            (.zmx ↔ .seq, by extension)
//! tpt-optics schema                        (JSON Schema of the design file)
//! tpt-optics mcp                           (JSON-RPC over stdio)
//! ```

mod convert;
mod design;
mod manifest;
mod mcp;

use std::fmt::Write as _;
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use serde_json::json;
use tpt_optics::core::index::IndexAt;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::geometric::sequential::field_ray;

use crate::design::{stop_text, Design};

const USAGE: &str =
    "usage: tpt-optics <info|glass|analyze|optimize|trace|run|convert|schema|mcp> ...";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(output) => {
            print!("{output}");
            ExitCode::SUCCESS
        }
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

/// Run one command and return what it prints.
fn run(args: &[String]) -> Result<String, String> {
    let (command, rest) = args.split_first().ok_or_else(|| USAGE.to_string())?;
    match command.as_str() {
        "info" => Ok(info()),
        "glass" => glass(rest),
        "analyze" => analyze(rest),
        "optimize" => optimize(rest),
        "trace" => trace(rest),
        "run" => run_design(rest),
        "convert" => convert_file(rest),
        "schema" => schema(),
        "mcp" => mcp_loop().map(|()| String::new()),
        other => Err(format!("unknown command '{other}'\n{USAGE}")),
    }
}

/// JSON Schema of the design file, for editors and validators.
fn schema() -> Result<String, String> {
    let schema = schemars::schema_for!(Design);
    serde_json::to_string_pretty(&schema)
        .map(|s| format!("{s}\n"))
        .map_err(|e| e.to_string())
}

fn info() -> String {
    format!(
        "tpt-optics {}\ncrates: core, geometric, wave, polarization, waveguide, design, io, verify\n\
         commands: info glass analyze optimize trace run convert mcp\n",
        env!("CARGO_PKG_VERSION")
    )
}

/// Value of `--flag VALUE` in `args`, if present.
fn flag_value<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

fn positional(args: &[String]) -> Option<&str> {
    args.iter()
        .find(|a| !a.starts_with("--"))
        .map(String::as_str)
}

fn read_design(path: &str) -> Result<(Design, Vec<u8>), String> {
    let bytes = fs::read(path).map_err(|e| format!("cannot read '{path}': {e}"))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| format!("'{path}' is not UTF-8: {e}"))?;
    let design = Design::from_toml(text).map_err(|e| format!("'{path}': {e}"))?;
    Ok((design, bytes))
}

fn glass(args: &[String]) -> Result<String, String> {
    let name = positional(args).ok_or("usage: tpt-optics glass <NAME> [--nm <WAVELENGTH>]")?;
    let nm: f64 = flag_value(args, "--nm")
        .map(|v| v.parse().map_err(|_| format!("bad wavelength '{v}'")))
        .transpose()?
        .unwrap_or(587.6);
    let catalog = Design::catalog();
    let material = catalog
        .get(name)
        .ok_or_else(|| format!("glass '{name}' is not in the catalogue"))?;
    let wl = Wavelength::from_nm(nm);
    let n = material.index_at(wl, Kelvin::from_celsius(20.0)).n;
    let mut out = format!(
        "{}\n  n({nm} nm) = {n:.6}\n  Abbe number V_d = {:.2}\n  glass code = {}\n",
        material.name(),
        material.abbe_number(),
        material.glass_code()
    );
    if material.outside_validity(wl) {
        out.push_str("  warning: the wavelength is outside the fitted range of this glass\n");
    }
    Ok(out)
}

fn analyze(args: &[String]) -> Result<String, String> {
    let path = positional(args).ok_or("usage: tpt-optics analyze <DESIGN.toml>")?;
    let (design, _) = read_design(path)?;
    let report = design.analyze().map_err(|e| e.to_string())?;
    if has_flag(args, "--json") {
        return Ok(format!("{}\n", report.to_json()));
    }
    if has_flag(args, "--csv") {
        return Ok(report.to_csv());
    }
    if let Some(index) = flag_value(args, "--svg") {
        let i: usize = index
            .parse()
            .map_err(|_| format!("bad field index '{index}'"))?;
        return report
            .spot_svg(i)
            .ok_or_else(|| format!("no spot for field index {i}"));
    }
    let mut out = String::new();
    let fmt_opt =
        |v: Option<f64>| v.map_or_else(|| "n/a".to_string(), |x| format!("{:.4} mm", x * 1e3));
    let _ = writeln!(out, "design: {}", design.name);
    let _ = writeln!(out, "  EFL: {}", fmt_opt(report.efl));
    let _ = writeln!(out, "  BFL: {}", fmt_opt(report.bfl));
    let _ = writeln!(
        out,
        "  F/#: {}",
        report
            .f_number
            .map_or_else(|| "n/a".to_string(), |f| format!("{f:.2}"))
    );
    for f in &report.fields {
        let _ = writeln!(
            out,
            "  field {:>5.2}°: RMS spot {:.2} µm, distortion {:.4}, wavefront RMS {:.3e} m, best-focus shift {:.4} mm",
            f.angle_deg,
            f.rms_spot * 1e6,
            f.distortion,
            f.wavefront_rms,
            f.best_focus_shift * 1e3
        );
    }
    Ok(out)
}

fn optimize(args: &[String]) -> Result<String, String> {
    let path =
        positional(args).ok_or("usage: tpt-optics optimize <DESIGN.toml> [--write <OUT>]")?;
    let (design, _) = read_design(path)?;
    let (radii, result) = design.optimize().map_err(|e| e.to_string())?;
    let mut out = format!(
        "stop: {}\ncost: {:.6e}\nevaluations: {}\nradii (mm):\n",
        stop_text(result.stop),
        result.f,
        result.evaluations
    );
    for (surface, r) in design.surfaces.iter().zip(&radii) {
        let _ = writeln!(out, "  {r:.4}  (was {:.4})", surface.radius_mm);
    }
    if let Some(target) = flag_value(args, "--write") {
        let mut updated = design.clone();
        for (surface, r) in updated.surfaces.iter_mut().zip(&radii) {
            surface.radius_mm = *r;
        }
        let text = updated.to_toml().map_err(|e| e.to_string())?;
        fs::write(target, text).map_err(|e| format!("cannot write '{target}': {e}"))?;
        let _ = writeln!(out, "wrote {target}");
    }
    Ok(out)
}

fn trace(args: &[String]) -> Result<String, String> {
    let path = positional(args).ok_or("usage: tpt-optics trace <DESIGN.toml> --height-mm <H>")?;
    let (design, _) = read_design(path)?;
    let height: f64 = flag_value(args, "--height-mm")
        .ok_or("missing --height-mm")?
        .parse()
        .map_err(|_| "bad --height-mm".to_string())?;
    let field: f64 = flag_value(args, "--field-deg")
        .map(|v| v.parse().map_err(|_| format!("bad --field-deg '{v}'")))
        .transpose()?
        .unwrap_or(0.0);
    let sys = design.build().map_err(|e| e.to_string())?;
    let res = sys.trace(&field_ray(height * 1e-3, field), design.context());
    let mut out = String::new();
    for (i, p) in res.points.iter().enumerate() {
        let _ = writeln!(
            out,
            "surface {i}: x = {:.6} mm, y = {:.6} mm, z = {:.6} mm",
            p.x * 1e3,
            p.y * 1e3,
            p.z * 1e3
        );
    }
    if res.blocked {
        let _ = writeln!(out, "blocked: {:?}", res.block_reason);
    }
    let _ = writeln!(out, "OPD = {:.6e} m", res.opd);
    Ok(out)
}

fn run_design(args: &[String]) -> Result<String, String> {
    let path = positional(args).ok_or("usage: tpt-optics run <DESIGN.toml> --out <DIR>")?;
    let out_dir = PathBuf::from(flag_value(args, "--out").ok_or("missing --out <DIR>")?);
    let (design, bytes) = read_design(path)?;
    let key = manifest::run_key("run", &bytes);
    if let Some(cached) = manifest::load_cached(&out_dir, &key) {
        return Ok(format!("cache hit: {key}\n{cached}\n"));
    }
    let report = design.analyze().map_err(|e| e.to_string())?;
    let mut results = json!({ "analysis": serde_json::from_str::<serde_json::Value>(&report.to_json()).map_err(|e| e.to_string())? });
    if design.optimize.is_some() {
        let (radii, result) = design.optimize().map_err(|e| e.to_string())?;
        results["optimization"] = json!({
            "radii_mm": radii,
            "stop": stop_text(result.stop),
            "cost": result.f,
            "evaluations": result.evaluations,
            "history": result.history,
        });
    }
    let results_text = format!("{results:#}\n");
    write_outputs(&out_dir, &key, path, &results_text)?;
    manifest::store_cached(&out_dir, &key, &results_text)?;
    Ok(format!("cache miss: {key}\n{results_text}"))
}

fn write_outputs(out_dir: &Path, key: &str, input: &str, results: &str) -> Result<(), String> {
    fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;
    fs::write(out_dir.join("results.json"), results).map_err(|e| e.to_string())?;
    let manifest_text = format!(
        "{:#}\n",
        json!({
            "tool": "tpt-optics",
            "version": env!("CARGO_PKG_VERSION"),
            "command": "run",
            "input": input,
            "input_key": key,
            "hash_algorithm": "FNV-1a 64 (identifies inputs; not cryptographic)",
            "outputs": ["results.json"],
        })
    );
    fs::write(out_dir.join("manifest.json"), manifest_text).map_err(|e| e.to_string())
}

fn convert_file(args: &[String]) -> Result<String, String> {
    let input = args.first().ok_or("usage: tpt-optics convert <IN> <OUT>")?;
    let output = args.get(1).ok_or("usage: tpt-optics convert <IN> <OUT>")?;
    let text = fs::read_to_string(input).map_err(|e| format!("cannot read '{input}': {e}"))?;
    let converted = match (extension(input), extension(output)) {
        (Some("zmx"), Some("seq")) => convert::zmx_to_seq(&text)?,
        (Some("seq"), Some("zmx")) => convert::seq_to_zmx(&text)?,
        _ => {
            return Err("convert needs a .zmx or .seq input and the other format as output".into())
        }
    };
    fs::write(output, converted).map_err(|e| format!("cannot write '{output}': {e}"))?;
    Ok(format!("wrote {output}\n"))
}

fn extension(path: &str) -> Option<&str> {
    Path::new(path).extension().and_then(|e| e.to_str())
}

/// Serve JSON-RPC requests from stdin until it closes.
fn mcp_loop() -> Result<(), String> {
    let stdin = io::stdin();
    let mut stdout = io::stdout().lock();
    for line in stdin.lock().lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() {
            continue;
        }
        if let Some(response) = mcp::handle(&line) {
            writeln!(stdout, "{response}").map_err(|e| e.to_string())?;
            stdout.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
