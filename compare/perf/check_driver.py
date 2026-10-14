"""A check of the driver itself: that it watches the process that does the work.

    uv run python perf/check_driver.py

Starts, the way `run.py` starts the SciPy runner, a process that holds a
known amount of memory and goes through the phases of a runner, and holds
what the driver measured against it: the memory it read, and the core the
process was on. SciPy's memory was reported as 4 MB for every case from the
first comparison on, because the driver sampled a launcher and not the
interpreter (#166); nothing checked the driver.

Run before a comparison whose results are to be published, on the machine
that runs it.
"""

import json
import os
import sys
import time

HELD = 300 * 2**20

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--hold":
    # The process the driver starts: a runner that holds `HELD` bytes from
    # its second phase on.
    def mark(name: str) -> None:
        print(json.dumps({"mark": name, "t_ns": time.time_ns()}), flush=True)

    print(json.dumps({"mark": "pid", "pid": os.getpid()}), flush=True)
    import psutil

    time.sleep(0.2)
    mark("loaded")
    held = bytearray(HELD)
    held[:: 4096] = b"\x01" * len(held[:: 4096])
    time.sleep(0.2)
    mark("built")
    time.sleep(0.2)
    mark("evaluated")
    time.sleep(0.2)
    mark("timed")
    print(json.dumps({"mark": "result", "build_s": [0.1], "eval_s": [0.1], "cores": psutil.Process().cpu_affinity()}), flush=True)
    sys.exit(0)

import run  # noqa: E402  (next to this file)


def main() -> None:
    run.command = lambda library, directory: [sys.executable, "-X", "utf8", __file__, "--hold"]
    measured = run.run_one("scipy", run.HERE, pause=False)
    assert "failed" not in measured, measured
    rss = run.memory(measured)
    grown = rss["built"] - rss["loaded"]
    print(f"held {HELD / 2**20:.0f} MB; the driver read {rss['loaded'] / 2**20:.0f} MB before and {rss['built'] / 2**20:.0f} MB with it held")
    assert 0.9 * HELD <= grown <= 1.2 * HELD, f"the driver read {grown / 2**20:.0f} MB more, where {HELD / 2**20:.0f} MB was held"
    assert measured["cores"] == [run.CORE], f"the process was on cores {measured['cores']}, not on {run.CORE}"
    print(f"the process was on core {run.CORE}")
    print("ok")


if __name__ == "__main__":
    main()
