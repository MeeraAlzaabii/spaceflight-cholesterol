"""Download the processed GLDS-241 (OSD-241) RNA-seq tables from NASA OSDR.

- Files are saved unchanged in data/raw/OSD-241/ and made read-only.
- First run: writes metadata/checksums.sha256 and metadata/sources.tsv.
  Commit both. Later runs (and CI) only *check* the files against them, so a
  changed file on the server is caught instead of silently used.

Usage:
    python scripts/download_data.py           # download and check
    python scripts/download_data.py --list    # list what OSDR offers for the study
"""
import argparse
import hashlib
import json
import stat
import sys
import urllib.request
from datetime import date
from pathlib import Path

import yaml

BASE = "https://osdr.nasa.gov"


def find_files(obj):
    """Yield every dict in the OSDR API answer that has file_name and remote_url."""
    if isinstance(obj, dict):
        if "file_name" in obj and "remote_url" in obj:
            yield obj
        for value in obj.values():
            yield from find_files(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from find_files(value)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--list", action="store_true", help="list available files and exit")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))["data"]
    api = f"{BASE}/osdr/data/osd/files/{cfg['osd_id']}"
    raw_dir = Path(cfg["raw_dir"])
    checksums = Path(cfg["checksums"])
    sources = Path(cfg["sources"])

    with urllib.request.urlopen(api, timeout=120) as response:
        files = list(find_files(json.load(response)))
    if not files:
        sys.exit(f"No files found in the OSDR answer from {api}")

    if args.list:
        for f in sorted(files, key=lambda f: f["file_name"]):
            print(f["file_name"])
        return

    patterns = [p.lower() for p in cfg["download_patterns"]]
    wanted = sorted(
        (f for f in files
         if f["file_name"].lower().endswith(".csv")
         and any(p in f["file_name"].lower() for p in patterns)),
        key=lambda f: f["file_name"])
    if not wanted:
        sys.exit("No file matches download_patterns in config/config.yaml. "
                 "Run `pixi run python scripts/download_data.py --list` to see the names.")

    raw_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in wanted:
        url = f["remote_url"]
        if url.startswith("/"):
            url = BASE + url
        dest = raw_dir / f["file_name"]
        if dest.exists():
            print(f"already here: {dest}")
        else:
            print(f"downloading:  {f['file_name']}")
            part = dest.with_name(dest.name + ".part")
            urllib.request.urlretrieve(url, part)
            part.rename(dest)
            dest.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)  # read-only
        rows.append((dest, url))

    lines = [f"{sha256(dest)}  {dest.as_posix()}" for dest, _ in rows]

    if checksums.exists() and checksums.stat().st_size > 0:
        expected = {ln.split(None, 1)[1].strip(): ln.split(None, 1)[0]
                    for ln in checksums.read_text().splitlines() if ln.strip()}
        bad = [line.split("  ", 1)[1] for line in lines
               if expected.get(line.split("  ", 1)[1]) != line.split("  ", 1)[0]]
        if bad:
            sys.exit("Checksum mismatch (the files differ from the recorded ones):\n  "
                     + "\n  ".join(bad)
                     + f"\nIf the change is expected, delete {checksums} and "
                       f"{sources}, rerun, and commit the new files.")
        print(f"checksums match {checksums}")
    else:
        checksums.parent.mkdir(parents=True, exist_ok=True)
        checksums.write_text("\n".join(lines) + "\n")
        with open(sources, "w") as fh:
            fh.write("file\turl\tbytes\tdownloaded\n")
            for dest, url in rows:
                fh.write(f"{dest.as_posix()}\t{url}\t{dest.stat().st_size}\t"
                         f"{date.today().isoformat()}\n")
        print(f"FIRST RUN: wrote {checksums} and {sources}. Commit them to Git.")


if __name__ == "__main__":
    main()