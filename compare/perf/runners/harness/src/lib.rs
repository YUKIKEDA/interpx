//! What every Rust runner of the comparison shares.
//!
//! A runner is one process for one case and one library. It reads the case
//! from a directory the driver (`compare/perf/run.py`) wrote, and goes through
//! the same phases as every other runner, so that the driver can lay the
//! memory of the process, which it samples from outside, against them:
//!
//! 1. the data is read (`loaded`),
//! 2. the interpolator is built once and kept (`built`),
//! 3. the queries are evaluated once (`evaluated`), and the values are written
//!    for the driver to compare with the other libraries',
//! 4. construction and evaluation are each repeated and timed (`timed`).
//!
//! Each phase ends with a line of JSON on standard output that holds the time
//! of day, and the last line holds the times.

use serde_json::{json, Value};
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

/// A case: its parameters, and the directory its arrays are in.
pub struct Case {
    dir: PathBuf,
    /// The parameters (`case.json`).
    pub spec: Value,
}

impl Case {
    /// The case in the directory named by the first argument of the process.
    pub fn from_args() -> Self {
        let dir = PathBuf::from(std::env::args().nth(1).expect("usage: runner CASE_DIR"));
        let text = fs::read_to_string(dir.join("case.json")).expect("case.json is there");
        let spec = serde_json::from_str(&text).expect("case.json is JSON");
        Self { dir, spec }
    }

    /// The family of the case (`oned`, `grid`, ...).
    pub fn family(&self) -> &str {
        self.text("family")
    }

    /// A parameter that is a string.
    pub fn text(&self, key: &str) -> &str {
        self.spec[key].as_str().unwrap_or_else(|| panic!("{key} is a string"))
    }

    /// A parameter that is a count.
    pub fn count(&self, key: &str) -> usize {
        self.spec[key].as_u64().unwrap_or_else(|| panic!("{key} is a count")) as usize
    }

    /// A parameter that is a number.
    pub fn number(&self, key: &str) -> f64 {
        self.spec[key].as_f64().unwrap_or_else(|| panic!("{key} is a number"))
    }

    /// A parameter that is a count or absent.
    pub fn count_or_none(&self, key: &str) -> Option<usize> {
        self.spec[key].as_u64().map(|n| n as usize)
    }

    /// The array `name`: 64-bit floats in the byte order of the machine.
    pub fn array(&self, name: &str) -> Vec<f64> {
        read(&self.dir.join(format!("{name}.f64")))
    }

    /// The array `name` as rows of `D` numbers.
    pub fn rows<const D: usize>(&self, name: &str) -> Vec<[f64; D]> {
        let flat = self.array(name);
        assert_eq!(flat.len() % D, 0, "{name} holds rows of {D}");
        flat.chunks_exact(D).map(|row| core::array::from_fn(|d| row[d])).collect()
    }

    /// Where this library's values go.
    fn values_path(&self, library: &str) -> PathBuf {
        self.dir.join(format!("values.{library}.f64"))
    }
}

fn read(path: &Path) -> Vec<f64> {
    let bytes = fs::read(path).unwrap_or_else(|_| panic!("{} is there", path.display()));
    bytes
        .chunks_exact(8)
        .map(|chunk| f64::from_ne_bytes(chunk.try_into().expect("eight bytes")))
        .collect()
}

/// Says that the library has nothing for this case, and ends the process.
pub fn unsupported(why: &str) -> ! {
    println!("{}", json!({ "mark": "unsupported", "why": why }));
    std::process::exit(0);
}

fn mark(name: &str) {
    let now = SystemTime::now().duration_since(UNIX_EPOCH).expect("after 1970");
    println!("{}", json!({ "mark": name, "t_ns": now.as_nanos() as u64 }));
}

/// Long enough for the driver to sample the memory of a phase a few times.
const PAUSE: Duration = Duration::from_millis(60);
/// The time one of the two timed loops takes at least.
const BUDGET: Duration = Duration::from_millis(250);
const FEWEST: usize = 5;
const MOST: usize = 2000;

/// The shortest time that is measured as one sample. The clock of the machine
/// has steps of 100 ns (Windows), so a call shorter than this is repeated
/// inside a sample until the sample is this long.
const SAMPLE: Duration = Duration::from_micros(200);

/// The seconds one call of `run` takes, sampled several times: a warm-up
/// first, then samples until [`BUDGET`] is used, at least [`FEWEST`] of them
/// (three when one call alone takes longer than the budget). A sample is a
/// number of calls that takes [`SAMPLE`] or more, divided by that number.
fn time(mut run: impl FnMut()) -> Vec<f64> {
    let start = Instant::now();
    run();
    let first = start.elapsed();
    // How many calls make a sample, from a second look at a few of them.
    let mut calls = 1_u32;
    loop {
        let start = Instant::now();
        for _ in 0..calls {
            run();
        }
        if start.elapsed() >= SAMPLE || calls >= 1 << 24 {
            break;
        }
        calls *= 2;
    }
    let fewest = if first > BUDGET { 3 } else { FEWEST };
    let mut seconds = Vec::new();
    let begun = Instant::now();
    while seconds.len() < fewest || (begun.elapsed() < BUDGET && seconds.len() < MOST) {
        let start = Instant::now();
        for _ in 0..calls {
            run();
        }
        seconds.push(start.elapsed().as_secs_f64() / f64::from(calls));
    }
    seconds
}

/// Goes through the phases for one library: `build` makes the interpolator,
/// `eval` writes the value at every query into its second argument, which has
/// `values` numbers.
pub fn run<B>(
    case: &Case,
    library: &str,
    values: usize,
    build: impl Fn() -> B,
    eval: impl Fn(&B, &mut [f64]),
) {
    mark("loaded");
    std::thread::sleep(PAUSE);
    let built = build();
    mark("built");
    std::thread::sleep(PAUSE);
    let mut out = vec![0.0_f64; values];
    eval(&built, &mut out);
    mark("evaluated");
    std::thread::sleep(PAUSE);
    let bytes: Vec<u8> = out.iter().flat_map(|v| v.to_ne_bytes()).collect();
    fs::write(case.values_path(library), bytes).expect("the values are written");
    drop(built);

    let build_s = time(|| {
        black_box(build());
    });
    let built = build();
    let eval_s = time(|| {
        eval(black_box(&built), &mut out);
        black_box(&out);
    });
    mark("timed");
    println!("{}", json!({ "mark": "result", "build_s": build_s, "eval_s": eval_s }));
}
