"""The real-world data of the comparison (docs/design.ja.md section 8.4.7).

    uv run python perf/real.py            # fetch what is missing, and check it

Two data sets, each taken whole, in the unit its source gives it out in. What
is used of each, and how its points are divided into data and points held
out, was fixed before anything was measured (the design has the rules and why).
Nothing fetched is committed: it goes to `perf/out/real/`.

- Elevation: one tile of the 1 m elevation model of Auckland, of Toitū Te
  Whenua Land Information New Zealand (CC BY 4.0). A file that does not
  change, with the checksum its source publishes.
- Weather: the daily maximum temperature at every station of NOAA's
  GHCN-Daily on one day (CC0 1.0). Its files are written again every day and
  have no dated copies, so no checksum can be fixed: the version and the
  checksums of what was fetched are recorded, and a later fetch is told when
  they differ.
"""

import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REAL = HERE / "out" / "real"

# The tile of the 1:10 000 sheet layout that holds the summit of Maungawhau
# (Mount Eden), in the newest survey that has it. 4800 by 7200 cells of 1 m.
# The hill crosses the south edge of the tile; the tile is used as it is, and
# no region is cut around the hill.
ELEVATION = {
    "name": "auckland-ba32-0401",
    "file": "BA32_10000_0401.tiff",
    "url": "https://nz-elevation.s3.ap-southeast-2.amazonaws.com/auckland/auckland-part-1_2024/dem_1m/2193/BA32_10000_0401.tiff",
    # From the tile's STAC item (`file:checksum`, a multihash of SHA-256).
    "sha256": "06b2c8bd1cf34188e2a20886dfb1069083a31198d4ec7d2f8482cf1e94090b92",
    "title": "Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401",
    "credit": "Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).",
}

# January 1 of the last calendar year that is complete, and the element with
# the most stations after precipitation (which is zero at most of them).
WEATHER = {
    "name": "ghcn-tmax-20250101",
    "day": "20250101",
    "element": "TMAX",
    "base": "https://noaa-ghcn-pds.s3.amazonaws.com/",
    "files": {"2025.csv.gz": "csv.gz/by_year/2025.csv.gz", "ghcnd-stations.txt": "ghcnd-stations.txt", "ghcnd-version.txt": "ghcnd-version.txt"},
    "title": "GHCN-Daily, maximum temperature of 2025-01-01, every station",
    "credit": "NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).",
}

# Of every this many scattered points, in the order of the source, the last
# is held out.
HELD_OUT_EVERY = 10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    print(f"fetching {url}", flush=True)
    partial = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(url, partial)
    partial.replace(path)


def fetch() -> None:
    """Fetches what is missing, and checks what is there."""
    REAL.mkdir(parents=True, exist_ok=True)
    tile = REAL / ELEVATION["file"]
    if not tile.exists():
        download(ELEVATION["url"], tile)
    found = sha256(tile)
    if found != ELEVATION["sha256"]:
        sys.exit(f"{tile.name}: SHA-256 {found}, where the source publishes {ELEVATION['sha256']}. Delete the file and fetch again.")
    print(f"{tile.name}: {tile.stat().st_size} bytes, checksum as published")

    fetched = False
    for name, key in WEATHER["files"].items():
        if not (REAL / name).exists():
            download(WEATHER["base"] + key, REAL / name)
            fetched = True
    record = REAL / "weather.json"
    now = {
        "version": (REAL / "ghcnd-version.txt").read_text(encoding="utf-8").splitlines()[0],
        "sha256": {name: sha256(REAL / name) for name in WEATHER["files"]},
    }
    if record.exists() and not fetched:
        before = json.loads(record.read_text(encoding="utf-8"))
        if before != now:
            sys.exit("the weather files are not the ones recorded in weather.json: delete perf/out/real/ and fetch again")
    record.write_text(json.dumps(now, indent=1), encoding="utf-8")
    print(f"weather: {now['version']}")
    for name, digest in now["sha256"].items():
        print(f"  {name}: {digest}")


def available() -> bool:
    return (REAL / ELEVATION["file"]).exists() and (REAL / "weather.json").exists()


def elevation() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The northings (ascending), the eastings, and the elevation at each pair
    of them, in metres: the centres of the cells of the tile."""
    import rasterio

    with rasterio.open(REAL / ELEVATION["file"]) as source:
        band = source.read(1).astype(np.float64)
        transform = source.transform
        assert source.nodata is None or not np.any(band == source.nodata), "the tile has no cell without a value"
    assert transform.a == 1.0 and transform.e == -1.0 and transform.b == 0.0 and transform.d == 0.0, "cells of 1 m, north up"
    eastings = transform.c + 0.5 + np.arange(band.shape[1])
    northings = transform.f - 0.5 - np.arange(band.shape[0])
    # The rows of the file go south; the axis is given ascending.
    return northings[::-1].copy(), eastings, band[::-1].copy()


def kept(size: int) -> np.ndarray:
    """The places of the grid's own points along an axis of `size` cells:
    every second one, and the last, so that every cell left out lies inside
    the grid of those kept."""
    places = np.arange(0, size, 2)
    return places if places[-1] == size - 1 else np.append(places, size - 1)


def elevation_split() -> dict:
    """The tile divided: a grid of every second row and column (and the
    last ones), and every other cell as a point held out."""
    northings, eastings, band = elevation()
    rows, columns = kept(len(northings)), kept(len(eastings))
    data = np.zeros(band.shape, dtype=bool)
    data[np.ix_(rows, columns)] = True
    held_rows, held_columns = np.nonzero(~data)
    q = np.stack([northings[held_rows], eastings[held_columns]], axis=1)
    return {
        "axis0": northings[rows],
        "axis1": eastings[columns],
        "values": band[np.ix_(rows, columns)],
        "q": q,
        "truth": band[held_rows, held_columns],
    }


def weather() -> tuple[np.ndarray, np.ndarray, dict]:
    """The longitude and latitude (degrees) of every station with a value of
    the element on the day that passed the source's checks, and the values
    (degrees Celsius), in the order of the file. Of stations at the same
    coordinates the first is taken."""
    cache = REAL / f"{WEATHER['name']}.npz"
    if cache.exists():
        loaded = np.load(cache)
        return loaded["points"], loaded["values"], json.loads(str(loaded["counts"]))
    stations = {}
    for line in (REAL / "ghcnd-stations.txt").read_text(encoding="utf-8").splitlines():
        stations[line[0:11]] = (float(line[21:30]), float(line[12:20]))
    points, values, seen = [], [], set()
    counts = {"rows": 0, "failed_a_check": 0, "same_coordinates": 0}
    with gzip.open(REAL / next(iter(WEATHER["files"])), "rt", encoding="utf-8", newline="") as text:
        for station, day, element, value, _measurement, quality, *_ in csv.reader(text):
            if day != WEATHER["day"] or element != WEATHER["element"]:
                continue
            counts["rows"] += 1
            if quality:
                counts["failed_a_check"] += 1
                continue
            place = stations[station]
            if place in seen:
                counts["same_coordinates"] += 1
                continue
            seen.add(place)
            points.append(place)
            # Tenths of a degree.
            values.append(int(value) / 10.0)
    points, values = np.array(points), np.array(values)
    np.savez(cache, points=points, values=values, counts=json.dumps(counts))
    return points, values, counts


def weather_split() -> dict:
    """The stations divided: of every ten in the order of the file, the last
    is held out."""
    points, values, _ = weather()
    held = np.arange(len(points)) % HELD_OUT_EVERY == HELD_OUT_EVERY - 1
    return {"points": points[~held], "values": values[~held], "q": points[held], "truth": values[held]}


if __name__ == "__main__":
    fetch()
    grid = elevation_split()
    print(f"elevation: a grid of {grid['values'].shape[0]} by {grid['values'].shape[1]}, {len(grid['q'])} cells held out")
    points, values, counts = weather()
    scattered = weather_split()
    print(f"weather: {counts}, {len(points)} stations: {len(scattered['points'])} as data, {len(scattered['q'])} held out")
