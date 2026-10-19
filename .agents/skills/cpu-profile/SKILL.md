---
name: cpu-profile
description: >-
  Sample a Windows release binary with Visual Studio's CPU profiler and
  resolve stacks, including ntdll frames. Use when the user asks for a
  profile, a flame graph, a hot path, or why an interpx evaluation path is slow
  compared with SciPy or another Rust crate.
---

# CPU profile on this machine

The agent shell is medium integrity. Administrators is deny-only. A kernel logger does not start.

Do not run `wpr` or `xperf -on`. Do not switch to in-process timers because a session failed to start.

## Collect

Profile a release binary that runs an interpx workload (a `compare/perf` runner or a bench-like example) for long enough to sample. Build it with debug info so the PDB matches the image.

```powershell
$env:CARGO_PROFILE_RELEASE_DEBUG = "1"
cargo build --release
```

`VSDiagnostics.exe` is under Visual Studio's Diagnostics Hub collector. `CpuUsageHigh.json` is next to it in `AgentConfigs`. Session id is an integer from 1 to 255.

```powershell
$vs = "C:\Program Files\Microsoft Visual Studio\18\Community\Team Tools\DiagnosticsHub\Collector\VSDiagnostics.exe"
$cfg = "C:\Program Files\Microsoft Visual Studio\18\Community\Team Tools\DiagnosticsHub\Collector\AgentConfigs\CpuUsageHigh.json"
& $vs start 4 "/launch:$exe" "/launchArgs:$args" "/loadConfig:$cfg" "/scratchLocation:$scratch"
```

If the start says the class is not registered, start again with a different integer session id. Do not treat that as the profiler being unavailable.

When the process has exited:

```powershell
& $vs stop 4 "/output:$scratch\cpu.diagsession"
& $vs expandDiagSession "$scratch\cpu.diagsession"
```

The sample ETL is `sc.user_aux.etl` under the expanded directory. The agent samples at 4000 Hz.

## Decode

Point the symbol path at the directory that contains the PDB.

```powershell
$env:_NT_SYMBOL_PATH = $pdbDir
$env:_NT_SYMCACHE_PATH = "$scratch\symcache"
xperf -i $etl -symbols -o "$scratch\profile.txt" -a profile -detail
xperf -i $etl -symbols -o "$scratch\stacks.txt" -a stack -butterfly 50 -process $procName -pid $pid
```

`-o` comes before `-a`. `profile -detail` is the leaf hit. `stack -butterfly` is the inclusive tree. Read a function's own block, the one that contains `***itself***`. Rows above that block are the parent's other callees.

## ntdll

Frames shown as `ntdll.dll!***unknown***` still have an RVA. Resolve it to the nearest export at or below that RVA:

```powershell
dumpbin /exports C:\Windows\System32\ntdll.dll
```

`RtlAllocateHeap` and `RtlFreeHeap` are the usual hot exports. `RtlQueryPerformanceCounter` is `Instant::now`, which means a timer was compiled into the binary. Report that share separately from the evaluation path.

## Disassembly

To say whether a loop is a branch or a conditional move, read the instructions, not the profile. `dumpbin` is next to the MSVC compiler; call it from PowerShell (Git Bash turns `/nologo` into a path), on a binary whose PDB is beside it under the name cargo gave it:

```powershell
$dumpbin = Join-Path (Get-ChildItem "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC" | Select-Object -First 1).FullName "bin\Hostx64\x64\dumpbin.exe"
& $dumpbin /nologo /disasm:nobytes $exe /out:"$scratch\out.asm"
```

The names are Rust's v0 mangling (`_RINvNtCs..._7interpx6search15eval_batch_into...`); search them by the parts of the path. Split the file into functions at the lines that end in `:` and do not start with a space, and classify each short backward jump's loop by its instructions: a search by halves is a loop with `shr` and a compare of floats (`ucomisd`), and it has either a `cmov` or only jumps. A rebuild overwrites the PDB: copy the binary and its PDB together before building another version.

A guess about what the compiler did then (an Issue that says "it became a branch in one place") is checked with the compiler of then. The comparison records it (`rustc` in `compare/perf/out/environment*.json`). Install it beside the current one, build the code of then from a worktree, read, and remove both:

```text
rustup toolchain install 1.99.0 --profile minimal
cargo +1.99.0 build --release --manifest-path <worktree>/compare/perf/runners/interpx/Cargo.toml
rustup toolchain uninstall 1.99.0
```

In #194 the current compiler alone left "it may have been a branch then" open; 1.99.0 on the code before #137 showed 34 conditional moves and no branch.

## What to report

Inclusive shares of the process samples, for the functions that own the time, and the exclusive heap exports. A leaf list alone hides a caller that spends its time in callees. State the sample count, the rate, and the binary's commit.

## VTune

`vtune.exe` is under `C:\Program Files (x86)\Intel\oneAPI\vtune\latest\bin64`. Use it when the question is which stall, not only which function.

Profile the binary where cargo built it. A copy without its PDB shows `func@0x...` instead of names.

Pin the target, not the collector:

```powershell
& $vtune -collect hotspots -knob sampling-mode=sw -r $result -- cmd /c "start /affinity 10 /b /wait $exe $args"
```

User-mode sampling (`sampling-mode=sw`) runs from the agent shell. It gives shares by function and module at a 10 ms interval, so repeat the build until the run lasts about 10 s.

Hardware event-based sampling (`-collect uarch-exploration`, `-collect memory-access`) needs an elevated process. Write the collection as a script and start it with `Start-Process pwsh -Verb RunAs -Wait`; the owner approves the UAC prompt. Reports run unelevated.

```powershell
& $vtune -report summary -r $result
& $vtune -report hotspots -r $result -group-by function -format csv -csv-delimiter "|"
& $vtune -report gprof-cc -r $result -format text -report-output $file
& $vtune -report callstacks -r $result -filter "module=ntdll.dll" -format csv -csv-delimiter "|"
```

`summary` of a `uarch-exploration` result has the top-down split. `hotspots` by function has the same columns per function. `gprof-cc` has inclusive time. `callstacks` filtered to `ntdll.dll` or `VCRUNTIME140.dll` attributes heap and copy time to the first `interpx::` frame.

Inlining folds callees into their caller: a Horner loop and its interval search can both appear under `PPoly::eval`. Read line-level rows with that in mind.

To count allocations, give the scratch binary a `#[global_allocator]` that wraps `System` and counts `alloc` and `realloc`. The count per construction and per evaluated point says whether the heap share is many small allocations. An evaluation path should show zero (`.cursor/rules/perf.mdc`).

Delete the result directories with the scratch binary (`.cursor/rules/scratch.mdc`). State the finding on the Issue: where the time goes decides the order of speed rows (`.cursor/rules/bench.mdc`).
