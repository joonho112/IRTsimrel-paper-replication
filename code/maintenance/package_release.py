"""Create a source ZIP from the checked release manifest."""
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

from pathlib import Path
import csv
import hashlib
import zipfile

ROOT = Path(__file__).resolve().parents[2]
AUTHOR = b"JoonHo Lee (jlee296@ua.edu)"
manifest = ROOT / "manifest/files.csv"
with manifest.open(newline="", encoding="utf8") as stream:
    records = list(csv.DictReader(stream))
for record in records:
    path = ROOT / record["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise SystemExit("Changed release file: " + record["path"])
output = ROOT / "output/release/IRTsimrel-paper-replication.zip"
output.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
    archive.comment = b"Author: " + AUTHOR
    paths = [ROOT / record["path"] for record in records] + [manifest]
    for path in paths:
        info = zipfile.ZipInfo("IRTsimrel-paper-replication/" + path.relative_to(ROOT).as_posix())
        info.compress_type = zipfile.ZIP_DEFLATED
        info.comment = b"Author: " + AUTHOR
        info.external_attr = 0o100644 << 16
        archive.writestr(info, path.read_bytes(), compresslevel=9)
print(output.relative_to(ROOT))
print(hashlib.sha256(output.read_bytes()).hexdigest())
