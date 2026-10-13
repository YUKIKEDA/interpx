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
        self.spec[key]
            .as_str()
            .unwrap_or_else(|| panic!("{key} is a string"))
    }

    /// A parameter that is a count.
    pub fn count(&self, key: &str) -> usize {
        self.spec[key]
            .as_u64()
            .unwrap_or_else(|| panic!("{key} is a count")) as usize
    }

    /// A parameter that is a number.
    pub fn number(&self, key: &str) -> f64 {
        self.spec[key]
            .as_f64()
            .unwrap_or_else(|| panic!("{key} is a number"))
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
        flat.as_chunks::<D>().0.to_vec()
    }

    /// Where this library's values go.
    fn values_path(&self, library: &str) -> PathBuf {
        self.dir.join(format!("values.{library}.f64"))
    }
}

fn read(path: &Path) -> Vec<f64> {
    let bytes = fs::read(path).unwrap_or_else(|_| panic!("{} is there", path.display()));
    bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|chunk| f64::from_ne_bytes(*chunk))
        .collect()
}

/// Says that the library has nothing for this case, and ends the process.
pub fn unsupported(why: &str) -> ! {
    println!("{}", json!({ "mark": "unsupported", "why": why }));
    std::process::exit(0);
}

fn mark(name: &str) {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("after 1970");
    println!("{}", json!({ "mark": name, "t_ns": now.as_nanos() as u64 }));
}

/// The time one of the two timed loops takes at least. With samples of
/// [`SAMPLE`] this is 250 samples of a short call. (It was five times as
/// long at first, and a run of all the cases took two hours, most of it in
/// these loops and in the pauses.)
const BUDGET: Duration = Duration::from_millis(50);

/// The pause after each phase: long enough for the driver to sample the
/// memory of the phase a few times. The driver asks for it (`PERF_PAUSE_MS`)
/// for the cases whose memory it reports, and for no others.
fn pause() {
    let asked = std::env::var("PERF_PAUSE_MS")
        .ok()
        .and_then(|ms| ms.parse().ok());
    if let Some(ms) = asked {
        std::thread::sleep(Duration::from_millis(ms));
    }
}
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

/// [`time`] for a call that takes its input by value: the inputs of a sample
/// are made before the sample starts, so that the time is that of the call
/// alone.
fn time_given<I>(make: impl Fn() -> I, mut run: impl FnMut(I)) -> Vec<f64> {
    let start = Instant::now();
    run(make());
    let first = start.elapsed();
    let mut calls = 1_usize;
    loop {
        let inputs: Vec<I> = (0..calls).map(|_| make()).collect();
        let start = Instant::now();
        inputs.into_iter().for_each(&mut run);
        // The inputs of one sample are all in memory at once, so a sample
        // has no more than this many calls.
        if start.elapsed() >= SAMPLE || calls >= 1 << 10 {
            break;
        }
        calls *= 2;
    }
    let fewest = if first > BUDGET { 3 } else { FEWEST };
    let mut seconds = Vec::new();
    let begun = Instant::now();
    while seconds.len() < fewest || (begun.elapsed() < BUDGET && seconds.len() < MOST) {
        let inputs: Vec<I> = (0..calls).map(|_| make()).collect();
        let start = Instant::now();
        inputs.into_iter().for_each(&mut run);
        seconds.push(start.elapsed().as_secs_f64() / calls as f64);
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
    run_given(case, library, values, || (), |()| build(), eval);
}

/// [`run`] for a construction that takes its data by value: `make` gives the
/// data (a copy of what was read), outside the time of `build`.
pub fn run_given<I, B>(
    case: &Case,
    library: &str,
    values: usize,
    make: impl Fn() -> I,
    build: impl Fn(I) -> B,
    eval: impl Fn(&B, &mut [f64]),
) {
    mark("loaded");
    pause();
    let built = build(make());
    mark("built");
    pause();
    let mut out = vec![0.0_f64; values];
    eval(&built, &mut out);
    mark("evaluated");
    pause();
    let bytes: Vec<u8> = out.iter().flat_map(|v| v.to_ne_bytes()).collect();
    fs::write(case.values_path(library), bytes).expect("the values are written");
    drop(built);

    let build_s = time_given(&make, |input| {
        black_box(build(input));
    });
    let built = build(make());
    let eval_s = time(|| {
        eval(black_box(&built), &mut out);
        black_box(&out);
    });
    mark("timed");
    println!(
        "{}",
        json!({ "mark": "result", "build_s": build_s, "eval_s": eval_s })
    );
}
