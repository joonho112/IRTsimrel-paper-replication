# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

from pathlib import Path
import hashlib, json, math, csv, os

HERE = Path(__file__).resolve()
MS = HERE.parents[3]
CB4 = MS / "data" / "precomputed"
DERIVED = MS / "output" / "derived"
FIG = MS / "output" / "figures"
TAB = MS / "output" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SOURCES = {
    "T_s2_accuracy": CB4 / "tables/T_s2_accuracy.csv",
    "s2_calibration": CB4 / "results/s2/calibration.csv",
    "T_s3_m": CB4 / "tables/T_s3_m.csv",
    "s3_m_sensitivity": CB4 / "results/s3/m_sensitivity.csv",
    "T_s6_arm_a": CB4 / "tables/T_s6_arm_a.csv",
    "s6_arm_a_runs": CB4 / "results/s6/arm_a_runs.csv",
    "normal_calibration": CB4 / "software/normal_calibration.csv",
    "normal_item_form": CB4 / "software/normal_item_form.csv",
    "normal_roots": CB4 / "software/normal_roots.csv",
    "unreachable_case": CB4 / "software/unreachable_case.csv",
    "unreachable_warnings": CB4 / "software/unreachable_warnings.txt",
    "T_s4_correspondence": CB4 / "tables/T_s4_supported_correspondence.csv",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def provenance(name, inputs, outputs, extra=None):
    rec = {
        "author": "JoonHo Lee (jlee296@ua.edu)",
        "tool": name,
        "inputs": {k: {"path": str(p.relative_to(MS)), "sha256": sha256(p)} for k, p in inputs.items()},
        "outputs": [str(Path(o).relative_to(MS)) for o in outputs],
    }
    if extra:
        rec.update(extra)
    DERIVED.mkdir(parents=True, exist_ok=True)
    with open(DERIVED / f"{name}.provenance.json", "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=1, ensure_ascii=False)
    return rec


def _binom_cdf(k, n, p):
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(0, k + 1))


def _bisect(g, a, b, it=60):
    for _ in range(it):
        m = (a + b) / 2
        if g(a) * g(m) <= 0:
            b = m
        else:
            a = m
    return (a + b) / 2


def clopper_pearson(x, n, alpha=0.05):
    """Exact two-sided Clopper-Pearson interval (pure Python; no scipy dependency)."""
    lo = 0.0 if x == 0 else _bisect(lambda p: 1 - _binom_cdf(x - 1, n, p) - alpha / 2, 0, 1)
    hi = 1.0 if x == n else _bisect(lambda p: _binom_cdf(x, n, p) - alpha / 2, 0, 1)
    return lo, hi
