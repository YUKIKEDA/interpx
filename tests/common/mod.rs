//! Shared helpers for integration tests: the golden-file loader and the
//! tolerance comparison.
//!
//! Golden files are written by `compare/generate_*.py` and read here; the format
//! is described in `compare/common/golden.py`. `cargo test` never runs Python.

// Each integration-test crate compiles this module on its own and uses a subset
// of it. The allow keeps `-D warnings` quiet for the unused part
// (docs/design.ja.md section 8.1: one shared helper module for all golden tests).
#![allow(dead_code)]

use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer};
use std::path::PathBuf;

/// Provenance recorded by the generator.
#[derive(Debug, Deserialize)]
pub struct Meta {
    pub generator: String,
    pub scipy: String,
    pub numpy: String,
    pub seed: u64,
}

/// A golden file whose cases have the family's own shape `C`.
#[derive(Debug, Deserialize)]
pub struct GoldenFile<C> {
    pub meta: Meta,
    pub cases: Vec<C>,
}

/// An `f64` as written by the generator: a JSON number, or one of the strings
/// `"NaN"`, `"Infinity"`, `"-Infinity"`.
#[derive(Debug, Clone, Copy)]
pub struct Real(pub f64);

impl<'de> Deserialize<'de> for Real {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum Raw {
            Num(f64),
            Text(String),
        }
        match Raw::deserialize(d)? {
            Raw::Num(v) => Ok(Real(v)),
            Raw::Text(s) => match s.as_str() {
                "NaN" => Ok(Real(f64::NAN)),
                "Infinity" => Ok(Real(f64::INFINITY)),
                "-Infinity" => Ok(Real(f64::NEG_INFINITY)),
                other => Err(serde::de::Error::custom(format!(
                    "not a real number: {other:?}"
                ))),
            },
        }
    }
}

/// Unwraps a column of [`Real`] into plain `f64`.
pub fn reals(v: &[Real]) -> Vec<f64> {
    v.iter().map(|r| r.0).collect()
}

/// Loads `compare/goldens/<name>.json`.
///
/// # Errors
///
/// Returns a message naming the path when the file is missing or malformed.
/// A test cannot continue without its golden data, so it calls `.expect(..)`.
pub fn load<C: DeserializeOwned>(name: &str) -> Result<GoldenFile<C>, String> {
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "compare", "goldens"]
        .iter()
        .collect::<PathBuf>()
        .join(format!("{name}.json"));
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("cannot read golden file {}: {e}", path.display()))?;
    serde_json::from_str(&text)
        .map_err(|e| format!("cannot parse golden file {}: {e}", path.display()))
}

/// A fixed relative and absolute tolerance, written in the test that uses it
/// (`.cursor/rules/numeric-tests.mdc`). Never derive it from the constant the
/// test guards.
#[derive(Debug, Clone, Copy)]
pub struct Tol {
    pub rtol: f64,
    pub atol: f64,
}

/// The class for direct methods in f64 (docs/design.ja.md section 8.2).
pub const DIRECT_F64: Tol = Tol {
    rtol: 1e-12,
    atol: 1e-12,
};

/// Whether `got` matches `want` within `tol`: `|got - want| <= atol + rtol * |want|`.
///
/// Two NaNs match. A NaN and a number do not. Infinities match only when equal
/// in sign.
pub fn close(got: f64, want: f64, tol: Tol) -> bool {
    if got.is_nan() || want.is_nan() {
        return got.is_nan() && want.is_nan();
    }
    if got.is_infinite() || want.is_infinite() {
        return got.is_infinite()
            && want.is_infinite()
            && got.is_sign_positive() == want.is_sign_positive();
    }
    (got - want).abs() <= tol.atol + tol.rtol * want.abs()
}

/// Asserts that `got` matches `want` within `tol`, naming `what` on failure.
pub fn assert_close(got: f64, want: f64, tol: Tol, what: &str) {
    assert!(
        close(got, want, tol),
        "{what}: got {got:e}, want {want:e} (rtol {}, atol {})",
        tol.rtol,
        tol.atol
    );
}
