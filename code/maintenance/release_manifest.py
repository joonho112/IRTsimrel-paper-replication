"""Record SHA-256, size, and authorship for the distributed release files."""
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

from pathlib import Path
import csv
import hashlib

ROOT = Path(__file__).resolve().parents[2]
AUTHOR = "JoonHo Lee (jlee296@ua.edu)"
EXCLUDED = {".git", ".library", ".renv", ".quarto", "_book", "__pycache__", ".Rproj.user",
            "dev", "log", "logs", "private", "notes", "drafts", "tmp", "temp"}


def release_files():
    for path in sorted(ROOT.rglob("*")):
        rel = path.relative_to(ROOT)
        if not path.is_file() or any(part in EXCLUDED for part in rel.parts):
            continue
        if rel.parts[0] == "output" and rel.as_posix() != "output/README.md":
            continue
        if rel.parts[:2] == ("data", "external") and path.name != "README.md":
            continue
        if (rel.as_posix() == "manifest/files.csv"
                or path.name in {".DS_Store", ".Rhistory", ".RData", ".Ruserdata", ".env"}
                or path.name.startswith(".env.")
                or path.suffix in {".log", ".tmp", ".bak", ".swp", ".pyc"}):
            continue
        yield path


def main():
    rows = []
    for path in release_files():
        upstream = path.relative_to(ROOT).as_posix() in {
            "data/irw/irw_diff_pool.csv", "data/irw/LICENSE-MIT.txt"}
        rows.append({"path": path.relative_to(ROOT).as_posix(),
                     "bytes": path.stat().st_size,
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "author": AUTHOR,
                     "source": "https://github.com/itemresponsewarehouse/Rpkg/tree/cc96e459448ee9d268eb233881e949c1151cbdb6" if upstream else "",
                     "source_copyright": "Copyright (c) 2025 irw authors" if upstream else "",
                     "source_license": "MIT" if upstream else ""})
    with (ROOT / "manifest/files.csv").open("w", newline="", encoding="utf8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["path", "bytes", "sha256", "author", "source", "source_copyright", "source_license"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Manifest: {len(rows)} files")


if __name__ == "__main__":
    main()
