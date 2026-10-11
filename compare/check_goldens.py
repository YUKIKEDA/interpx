"""Check the generated goldens as a set. Run last by `just gen-goldens`.

Fails when the directory is over the packaging budget or when a file is left over
from a generator that no longer exists.
"""

import json
import sys

from common import golden


def main() -> None:
    here = golden.GOLDENS.parent
    for path in sorted(golden.GOLDENS.glob("*.json")):
        meta = json.loads(path.read_text(encoding="utf-8"))["meta"]
        if not (here / meta["generator"]).is_file():
            sys.exit(f"{path.name}: its generator {meta['generator']} does not exist")
    golden.check_total_size()


if __name__ == "__main__":
    main()
