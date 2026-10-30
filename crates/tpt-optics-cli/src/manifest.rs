//! Run manifests and content-addressed caching.
//!
//! A run is identified by a hash of its input bytes and its command, using
//! FNV-1a (64-bit). The hash is for identifying identical inputs and for cache
//! lookup; it is not a cryptographic digest. A cache entry is a file named
//! after the hash, so the same input and command always reuse one result.

use std::fs;
use std::path::{Path, PathBuf};

/// FNV-1a 64-bit hash of `bytes`.
#[must_use]
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Identity of a run: the hash of its command and input.
#[must_use]
pub fn run_key(command: &str, input: &[u8]) -> String {
    let mut bytes = command.as_bytes().to_vec();
    bytes.push(0);
    bytes.extend_from_slice(input);
    alloc_key(fnv1a(&bytes))
}

fn alloc_key(h: u64) -> String {
    format!("fnv1a64-{h:016x}")
}

/// Path of the cache entry for `key` under `out_dir`.
#[must_use]
pub fn cache_path(out_dir: &Path, key: &str) -> PathBuf {
    out_dir.join("cache").join(format!("{key}.json"))
}

/// Cached output for `key`, if one exists.
#[must_use]
pub fn load_cached(out_dir: &Path, key: &str) -> Option<String> {
    fs::read_to_string(cache_path(out_dir, key)).ok()
}

/// Store `output` as the cached result for `key`.
///
/// # Errors
/// Returns the I/O error text if the directory or file cannot be written.
pub fn store_cached(out_dir: &Path, key: &str, output: &str) -> Result<(), String> {
    let path = cache_path(out_dir, key);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::write(path, output).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv1a_matches_the_reference_vectors() {
        // Reference values of FNV-1a 64.
        assert_eq!(fnv1a(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a(b"a"), 0xaf63_dc4c_8601_ec8c);
    }

    #[test]
    fn the_same_command_and_input_give_the_same_key() {
        assert_eq!(run_key("run", b"x"), run_key("run", b"x"));
        assert_ne!(run_key("run", b"x"), run_key("analyze", b"x"));
        assert_ne!(run_key("run", b"x"), run_key("run", b"y"));
    }

    #[test]
    fn cache_round_trips_through_the_filesystem() {
        let dir = std::env::temp_dir().join(format!("tpt-cache-test-{}", std::process::id()));
        let key = run_key("run", b"design");
        assert_eq!(load_cached(&dir, &key), None);
        store_cached(&dir, &key, "{\"ok\":true}").unwrap();
        assert_eq!(load_cached(&dir, &key).as_deref(), Some("{\"ok\":true}"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
