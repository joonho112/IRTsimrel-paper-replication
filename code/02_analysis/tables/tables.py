# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

import statistics as st
from collections import defaultdict
from common import SOURCES, DERIVED, TAB, read_csv, write_csv, provenance, clopper_pearson

TAB.mkdir(exist_ok=True)
prov_inputs = {}

# ---------------------------------------------------------------- S2 rows
cal = read_csv(SOURCES["s2_calibration"]); prov_inputs["s2_calibration"] = SOURCES["s2_calibration"]
assert len(cal) == 748 and all(r["package_version"] == "0.3.0" for r in cal)
ok = [r for r in cal if r["eqc_status"] == "ok"]
assert len(ok) == 743 and all(int(r["eqc_root_count"]) == 1 for r in ok)
mae = st.mean(abs(float(r["eqc_delta"])) for r in ok); mx = max(abs(float(r["eqc_delta"])) for r in ok)
fail = [r for r in cal if r["eqc_status"] != "ok"]
assert len(fail) == 5 and all(r["eqc_root_status"] == "above_upper" for r in fail)

model_lab = {"rasch": "Rasch", "2pl": "2PL", "3pl_g20": "3PL, $g=.20$", "3pl_beta": "3PL, beta $g$"}
shape_lab = {"normal": "normal", "bimodal": "bimodal", "skew_pos": "skewed", "heavy_tail": "heavy-tailed"}
pool_lab = {"parametric": "parametric", "irw": "IRW"}

# c* by shape at target .70, 20 items
grid = {}
for r in ok:
    if r["target_rho"] == "0.7" and r["n_items"] == "20":
        grid[(r["model"], r["pool"], r["latent_shape"])] = float(r["eqc_c_star"])
assert len(grid) == 32
lines = [r"\begin{tabular}{llrrrr}", r"\hline",
         r"Item model & Difficulty pool & normal & bimodal & skewed & heavy-tailed \\", r"\hline"]
csv_rows = []
for m in ["rasch", "2pl", "3pl_g20", "3pl_beta"]:
    for p in ["parametric", "irw"]:
        vals = [grid[(m, p, s)] for s in ["normal", "bimodal", "skew_pos", "heavy_tail"]]
        lines.append(f"{model_lab[m]} & {pool_lab[p]} & " + " & ".join(f"{v:.3f}" for v in vals) + r" \\")
        csv_rows.append({"model": m, "pool": p, **{s: grid[(m, p, s)] for s in ["normal", "bimodal", "skew_pos", "heavy_tail"]}})
lines += [r"\hline", r"\end{tabular}"]
(TAB / "tab_sc03_cstar_by_shape.tex").write_text("\n".join(lines) + "\n")
write_csv(DERIVED / "tab_sc03_cstar_by_shape.csv", csv_rows)

# failures
lines = [r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{llllrrl}", r"\hline",
         r"Item model & Latent shape & Pool & $I$ & Target & Largest $\tilde\rho$ found & Status \\", r"\hline"]
frows = []
for r in sorted(fail, key=lambda r: (r["model"], r["latent_shape"], r["pool"], int(r["n_items"]))):
    lines.append(f"{model_lab[r['model']]} & {shape_lab[r['latent_shape']]} & {pool_lab[r['pool']]} & {r['n_items']} & {float(r['target_rho']):.2f} & {float(r['eqc_rho_max']):.3f} & above $c_U$ \\\\")
    frows.append({k: r[k] for k in ("cond_id", "model", "latent_shape", "pool", "n_items", "target_rho", "eqc_rho_max", "eqc_root_status")})
lines += [r"\hline", r"\end{tabular}"]
(TAB / "tab_sc03_failures.tex").write_text("\n".join(lines) + "\n")
write_csv(DERIVED / "tab_sc03_failures.csv", frows)

# three-index gaps at the EQC root (used in F04 / Supplement A text)
gp = sorted(float(r["rho_tilde_at_cstar"]) - float(r["rho_psd_at_cstar"]) for r in ok if r["rho_psd_at_cstar"] not in ("NA", ""))
gw = sorted(float(r["rho_tilde_at_cstar"]) - float(r["w_bar_at_cstar"]) for r in ok if r["w_bar_at_cstar"] not in ("NA", ""))
gaps = {"psd_n": len(gp), "psd_median": st.median(gp), "psd_p90": gp[int(.9 * len(gp))], "psd_max": gp[-1],
        "wbar_n": len(gw), "wbar_median": st.median(gw), "wbar_p90": gw[int(.9 * len(gw))], "wbar_max": gw[-1]}
write_csv(DERIVED / "s2_gaps_at_cstar.csv", [gaps])

# ---------------------------------------------------------------- S3 holdout by M
ms = [r for r in read_csv(SOURCES["s3_m_sensitivity"]) if r["status"] == "ok"]; prov_inputs["s3_m_sensitivity"] = SOURCES["s3_m_sensitivity"]
byM = defaultdict(list)
for r in ms:
    byM[int(r["M"])].append(r)
lines = [r"\begin{tabular}{rrrrrr}", r"\hline",
         r"$M$ & runs & solver residual MAE & fresh-draw MAE & SD of signed fresh error & SD of $c^*$ \\", r"\hline"]
for M in sorted(byM):
    rr = byM[M]
    lines.append(f"{M:,} & {len(rr)} & {st.mean(abs(float(r['calibration_residual'])) for r in rr):.1e} & "
                 f"{st.mean(abs(float(r['holdout_error'])) for r in rr):.4f} & {st.pstdev([float(r['holdout_error']) for r in rr]):.4f} & "
                 f"{st.stdev(float(r['c_star']) for r in rr):.3f} \\\\")
lines += [r"\hline", r"\end{tabular}"]
(TAB / "tab_sc04_holdout_by_M.tex").write_text("\n".join(lines).replace("e-0", r"$\times 10^{-") .replace("e-1", r"$\times 10^{-1") + "\n")
# (the replace above is cosmetic only for exponents like 2.4e-09 -> handled below more robustly)
def fmt_e(x):
    m, e = f"{x:.1e}".split("e"); return f"${m}\\times10^{{{int(e)}}}$"
lines = [r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{rrrrrr}", r"\hline",
         r"$M$ & runs & residual MAE & fresh-draw MAE & SD of signed fresh error & SD of $c^*$ \\", r"\hline"]
for M in sorted(byM):
    rr = byM[M]
    lines.append(f"{M:,} & {len(rr)} & {fmt_e(st.mean(abs(float(r['calibration_residual'])) for r in rr))} & "
                 f"{st.mean(abs(float(r['holdout_error'])) for r in rr):.4f} & {st.pstdev([float(r['holdout_error']) for r in rr]):.4f} & "
                 f"{st.stdev(float(r['c_star']) for r in rr):.3f} \\\\")
lines += [r"\hline", r"\end{tabular}"]
(TAB / "tab_sc04_holdout_by_M.tex").write_text("\n".join(lines) + "\n")

# ---------------------------------------------------------------- Arm A tables
T = read_csv(SOURCES["T_s6_arm_a"]); prov_inputs["T_s6_arm_a"] = SOURCES["T_s6_arm_a"]
targets = ["0.5", "0.6", "0.7", "0.8", "0.887", "0.9"]
conds = [("0", "0"), ("0", "0.2"), ("0", "0.4"), ("0.4", "0"), ("0.4", "0.2"), ("0.4", "0.4")]
cond_lab = [r"mean 0, SD 0", r"mean 0, SD .2", r"mean 0, SD .4", r"mean .4, SD 0", r"mean .4, SD .2", r"mean .4, SD .4"]
def arm_table(scale, method, fname, with_ci):
    rows = {(r["target_rho"], r["te_mean"], r["te_sd"]): r for r in T if r["te_scale"] == scale and r["method"] == method}
    nreps = {int(r["n_reps"]) for r in rows.values()}; assert len(nreps) == 1; n = nreps.pop()
    lines = [r"\setlength{\tabcolsep}{3.5pt}", r"\begin{tabular}{l" + "r" * 6 + "}", r"\hline",
             r"& \multicolumn{3}{c}{reference effect mean 0} & \multicolumn{3}{c}{reference effect mean .4} \\",
             r"Target & SD 0 & SD .2 & SD .4 & SD 0 & SD .2 & SD .4 \\", r"\hline"]
    csv_rows = []
    for t in targets:
        cells = []
        for (m, s) in conds:
            r = rows[(t, m, s)]; p = float(r["reject_rate"]); x = round(p * n)
            if with_ci and ((m, s) in (("0.4", "0"), ("0", "0.4"))):
                lo, hi = clopper_pearson(x, n); cells.append(f"{p:.3f} [{lo:.3f}, {hi:.3f}]".replace("0.", "."))
            else:
                cells.append(f"{p:.3f}".replace("0.", "."))
            csv_rows.append({"te_scale": scale, "method": method, "target_rho": t, "te_mean": m, "te_sd": s,
                             "n_reps": n, "reject_rate": p, "se_ratio": r["se_ratio"], "mean_estimate": r["mean_estimate"]})
        lines.append(f"{float(t):.3f}".rstrip("0").replace("0.", ".") + " & " + " & ".join(cells) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    (TAB / fname).write_text("\n".join(lines) + "\n")
    write_csv(DERIVED / fname.replace(".tex", ".csv"), csv_rows)
    return n
n_theta = arm_table("theta_fixed", "Constant", "tab_se02_arm_a_theta.tex", True)
n_vary = arm_table("theta_fixed", "Varying", "tab_se02_arm_a_varying.tex", False)
n_logit = arm_table("logit_fixed", "Constant", "tab_se03_arm_a_logit.tex", False)
assert n_theta == 1000 and n_vary == 1000 and n_logit == 100
# calibrated common slope by target (from raw runs)
runs = read_csv(SOURCES["s6_arm_a_runs"]); prov_inputs["s6_arm_a_runs"] = SOURCES["s6_arm_a_runs"]
disc = defaultdict(list)
for r in runs:
    if r["te_scale"] == "theta_fixed" and r["method"] == "Constant" and r["te_mean"] == "0.4" and r["te_sd"] == "0":
        disc[r["target_rho"]].append(float(r["disc"]))
disc_rows = [{"target_rho": t, "n": len(v), "disc_mean": st.mean(v), "disc_min": min(v), "disc_max": max(v),
              "implied_logit_effect": 0.4 * st.mean(v) / 1.7} for t, v in sorted(disc.items(), key=lambda kv: float(kv[0]))]
write_csv(DERIVED / "arm_a_common_slope_by_target.csv", disc_rows)

summary = {"s2": {"n_cells": len(cal), "n_ok": len(ok), "mae": mae, "max_abs": mx, "n_fail": len(fail), "gaps": gaps},
           "s3": {str(M): {"n": len(byM[M]), "holdout_mae": st.mean(abs(float(r["holdout_error"])) for r in byM[M]),
                           "holdout_sd": st.pstdev([float(r["holdout_error"]) for r in byM[M]])} for M in sorted(byM)},
           "arm_a_disc": disc_rows}
provenance("tables", prov_inputs, sorted(str(p) for p in TAB.glob("*.tex")), {"summary": summary})
print("S2:", len(ok), f"{mae:.3e}", f"{mx:.3e}", "fail", len(fail)); print("gaps:", {k: round(v, 4) if isinstance(v, float) else v for k, v in gaps.items()})
print("disc:", [(d["target_rho"], round(d["disc_mean"], 3), round(d["implied_logit_effect"], 3)) for d in disc_rows])
