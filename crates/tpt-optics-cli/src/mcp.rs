//! A small JSON-RPC 2.0 server over stdio, speaking the Model Context
//! Protocol's tool interface.
//!
//! Each input line is one request, and each response is one line. Requests
//! without an `id` are notifications and get no response. Tools:
//!
//! - `glass_index` — `{name, wavelength_nm?}`: index, Abbe number, glass code,
//!   and whether the wavelength is outside the glass's fitted range.
//! - `analyze` — `{design_toml}`: the analysis report as JSON.
//! - `optimize` — `{design_toml}`: the optimised radii (mm) and stop reason.

use serde_json::{json, Value};
use tpt_optics::core::index::IndexAt;
use tpt_optics::core::{Kelvin, Wavelength};

use crate::design::{stop_text, Design};

const PROTOCOL_VERSION: &str = "2024-11-05";

/// Handle one request line. Returns `None` for notifications.
#[must_use]
pub fn handle(line: &str) -> Option<String> {
    let request: Value = match serde_json::from_str(line) {
        Ok(v) => v,
        Err(e) => {
            return Some(error_response(
                &Value::Null,
                -32700,
                &format!("parse error: {e}"),
            ))
        }
    };
    let id = request.get("id").cloned()?;
    let method = request.get("method").and_then(Value::as_str).unwrap_or("");
    let params = request.get("params").cloned().unwrap_or(Value::Null);
    let outcome = match method {
        "initialize" => Ok(json!({
            "protocolVersion": PROTOCOL_VERSION,
            "capabilities": { "tools": {} },
            "serverInfo": { "name": "tpt-optics", "version": env!("CARGO_PKG_VERSION") },
        })),
        "tools/list" => Ok(tools_list()),
        "tools/call" => Ok(call_tool(&params)),
        "ping" => Ok(json!({})),
        other => Err((-32601, format!("method not found: {other}"))),
    };
    Some(match outcome {
        Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }).to_string(),
        Err((code, message)) => error_response(&id, code, &message),
    })
}

fn error_response(id: &Value, code: i64, message: &str) -> String {
    json!({
        "jsonrpc": "2.0",
        "id": id,
        "error": { "code": code, "message": message },
    })
    .to_string()
}

fn tools_list() -> Value {
    json!({
        "tools": [
            {
                "name": "glass_index",
                "description": "Refractive index, Abbe number and glass code of a catalogue glass.",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "name": { "type": "string" },
                        "wavelength_nm": { "type": "number" }
                    },
                    "required": ["name"]
                }
            },
            {
                "name": "analyze",
                "description": "Analyse a lens design given as TOML text; returns a JSON report.",
                "inputSchema": {
                    "type": "object",
                    "properties": { "design_toml": { "type": "string" } },
                    "required": ["design_toml"]
                }
            },
            {
                "name": "optimize",
                "description": "Optimise the variables of a TOML design; returns the radii (mm).",
                "inputSchema": {
                    "type": "object",
                    "properties": { "design_toml": { "type": "string" } },
                    "required": ["design_toml"]
                }
            }
        ]
    })
}

/// Result of a tool call: a text content block, with `isError` on failure.
fn tool_result(text: impl Into<String>, is_error: bool) -> Value {
    json!({
        "content": [{ "type": "text", "text": text.into() }],
        "isError": is_error,
    })
}

fn call_tool(params: &Value) -> Value {
    let name = params.get("name").and_then(Value::as_str).unwrap_or("");
    let args = params.get("arguments").cloned().unwrap_or(Value::Null);
    match name {
        "glass_index" => glass_index(&args),
        "analyze" => design_tool(&args, |d| d.analyze().map(|r| r.to_json())),
        "optimize" => design_tool(&args, |d| {
            d.optimize().map(|(radii, result)| {
                json!({
                    "radii_mm": radii,
                    "stop": stop_text(result.stop),
                    "cost": result.f,
                    "evaluations": result.evaluations,
                })
                .to_string()
            })
        }),
        other => tool_result(format!("unknown tool '{other}'"), true),
    }
}

fn glass_index(args: &Value) -> Value {
    let Some(name) = args.get("name").and_then(Value::as_str) else {
        return tool_result(String::from("missing string argument 'name'"), true);
    };
    let nm = args
        .get("wavelength_nm")
        .and_then(Value::as_f64)
        .unwrap_or(587.6);
    let catalog = Design::catalog();
    let Some(material) = catalog.get(name) else {
        return tool_result(format!("glass '{name}' is not in the catalogue"), true);
    };
    let wl = Wavelength::from_nm(nm);
    let n = material.index_at(wl, Kelvin::from_celsius(20.0)).n;
    let text = json!({
        "name": material.name(),
        "wavelength_nm": nm,
        "n": n,
        "abbe_number": material.abbe_number(),
        "glass_code": material.glass_code(),
        "outside_validity": material.outside_validity(wl),
    })
    .to_string();
    tool_result(text, false)
}

fn design_tool<F>(args: &Value, run: F) -> Value
where
    F: FnOnce(&Design) -> Result<String, crate::design::DesignError>,
{
    let Some(text) = args.get("design_toml").and_then(Value::as_str) else {
        return tool_result(String::from("missing string argument 'design_toml'"), true);
    };
    match Design::from_toml(text).and_then(|d| run(&d)) {
        Ok(out) => tool_result(out, false),
        Err(e) => tool_result(e.to_string(), true),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SINGLET: &str = r#"
[[surfaces]]
radius_mm = 100.0
semi_diameter_mm = 12.5
thickness_mm = 5.0
glass = "N-BK7"

[[surfaces]]
radius_mm = -100.0
semi_diameter_mm = 12.5
thickness_mm = 95.0

[[surfaces]]
radius_mm = 0.0
semi_diameter_mm = 12.5
thickness_mm = 0.0
"#;

    fn call(name: &str, arguments: &Value) -> String {
        let req = json!({
            "jsonrpc": "2.0", "id": 7, "method": "tools/call",
            "params": { "name": name, "arguments": arguments }
        });
        handle(&req.to_string()).unwrap()
    }

    fn text_of(response: &str) -> (Value, bool) {
        let v: Value = serde_json::from_str(response).unwrap();
        let result = &v["result"];
        let text = result["content"][0]["text"].as_str().unwrap().to_string();
        let is_error = result["isError"].as_bool().unwrap();
        (
            serde_json::from_str(&text).unwrap_or(Value::String(text)),
            is_error,
        )
    }

    #[test]
    fn initialize_reports_the_protocol_and_tools_capability() {
        let r = handle(r#"{"jsonrpc":"2.0","id":1,"method":"initialize"}"#).unwrap();
        let v: Value = serde_json::from_str(&r).unwrap();
        assert_eq!(v["result"]["protocolVersion"], PROTOCOL_VERSION);
        assert!(v["result"]["capabilities"]["tools"].is_object());
    }

    #[test]
    fn notifications_get_no_response() {
        assert!(handle(r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#).is_none());
    }

    #[test]
    fn unknown_methods_and_bad_json_are_errors() {
        let r: Value =
            serde_json::from_str(&handle(r#"{"jsonrpc":"2.0","id":2,"method":"nope"}"#).unwrap())
                .unwrap();
        assert_eq!(r["error"]["code"], -32601);
        let bad: Value = serde_json::from_str(&handle("{not json").unwrap()).unwrap();
        assert_eq!(bad["error"]["code"], -32700);
    }

    #[test]
    fn tools_list_names_the_three_tools() {
        let r = handle(r#"{"jsonrpc":"2.0","id":3,"method":"tools/list"}"#).unwrap();
        let v: Value = serde_json::from_str(&r).unwrap();
        let names: Vec<&str> = v["result"]["tools"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| t["name"].as_str().unwrap())
            .collect();
        assert_eq!(names, ["glass_index", "analyze", "optimize"]);
    }

    #[test]
    fn glass_index_returns_bk7_properties() {
        let (v, is_error) = text_of(&call("glass_index", &json!({"name": "N-BK7"})));
        assert!(!is_error);
        assert!((v["n"].as_f64().unwrap() - 1.5168).abs() < 1e-3, "{v}");
        assert_eq!(v["outside_validity"], false);
    }

    #[test]
    fn unknown_glass_is_a_tool_error() {
        let (text, is_error) = text_of(&call("glass_index", &json!({"name": "ZZZ"})));
        assert!(is_error);
        assert!(text.to_string().contains("ZZZ"));
    }

    #[test]
    fn analyze_returns_the_report() {
        let (v, is_error) = text_of(&call("analyze", &json!({"design_toml": SINGLET})));
        assert!(!is_error, "{v}");
        assert!(v["efl"].as_f64().unwrap() > 0.0);
        assert!(v["fields"].is_array());
    }
}
