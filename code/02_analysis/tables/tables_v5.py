# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

import statistics as st
from collections import defaultdict
from common import CB4, TAB, DERIVED, SOURCES, read_csv, write_csv, provenance, clopper_pearson

EXTRA = {
    "turning_points": CB4 / "results/theory/turning_points.csv",
    "guttman_limit": CB4 / "results/theory/guttman_limit.csv",
    "tail_truncation": CB4 / "results/theory/tail_truncation_sweep.csv",
    "T_s2_runtime": CB4 / "tables/T_s2_runtime.csv",
}
TAB.mkdir(exist_ok=True)

def fmt(x, d=3):
    s = f"{x:.{d}f}"
    return s.replace("0.", ".", 1) if s.startswith("0.") or s.startswith("-0.") else s

def sci(x):
    m, e = f"{x:.2e}".split("e"); return f"${m}\\times10^{{{int(e)}}}$"

shape_lab = {"normal": "normal", "skew_pos": "skewed", "bimodal": "bimodal", "heavy_tail": "heavy-tailed"}

# ---------------------------------------------------------------- turning points
tp = read_csv(EXTRA["turning_points"]); assert len(tp) == 32
assert all(r["tilde_mono_wide"] == "TRUE" for r in tp)
lines = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{Turning points of the pointwise and inverse-information curves on the 32 finite grids of \cref{fig:functionals}}",
         r"\label{tab:turning}", r"\small", r"\setlength{\tabcolsep}{4pt}",
         r"\begin{tabular}{llrrrrr}", r"\toprule",
         r" & & \multicolumn{2}{c}{pointwise $\rhopw$} & \multicolumn{2}{c}{inverse-information $\wbar$} & \\",
         r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
         r"Shape & $I$ & $c$ at maximum & maximum & $c$ at maximum & maximum & $\rhotilde$ monotone \\", r"\midrule"]
rows_csv = []
for g in ("0", "0.2"):
    lines.append(r"\multicolumn{7}{l}{\emph{" + ("2PL" if g == "0" else "3PL, $g=.20$") + r"}} \\")
    for sh in ("normal", "skew_pos", "bimodal", "heavy_tail"):
        for I in ("10", "20", "40", "60"):
            r = [x for x in tp if x["shape"] == sh and x["n_items"] == I and x["guessing"] == g][0]
            wb = "---" if sh == "heavy_tail" else fmt(float(r["wbar_max"]))
            wc = "---" if sh == "heavy_tail" else f"{float(r['c_turn_wbar']):.2f}"
            lines.append(f"{shape_lab[sh]} & {I} & {float(r['c_turn_psd']):.2f} & {fmt(float(r['psd_max']))} & {wc} & {wb} & yes \\\\")
            rows_csv.append({k: r[k] for k in r})
lines += [r"\bottomrule", r"\end{tabular}", "", r"\smallskip", r"\begin{minipage}{0.95\textwidth}", r"\footnotesize",
          r"\emph{Note.} Each grid evaluates the three indices on 4,001 density-weighted nodes with represented variance one at 91 log-spaced multipliers from .02 to 50. The maximum and its location are grid values, not population quantities. For the heavy-tailed population the finite-grid inverse-information curve is a truncation of a divergent expectation (\cref{prop:tail}) and is not reported. The mean-information curve increased across every grid.",
          r"\end{minipage}", r"\end{table}"]
(TAB / "tab_turning_points.tex").write_text("\n".join(lines) + "\n")
write_csv(DERIVED / "tab_turning_points.csv", rows_csv)
tp_psd = [float(r["c_turn_psd"]) for r in tp]; tp_w = [float(r["c_turn_wbar"]) for r in tp]

# ---------------------------------------------------------------- Guttman limit
gl = read_csv(EXTRA["guttman_limit"]); assert len(gl) >= 7
lines = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{One 10-item 2PL form under a standard normal population at increasing multipliers}",
         r"\label{tab:guttman}", r"\small", r"\begin{tabular}{rrrrrr}", r"\toprule",
         r"$c$ & mean information $A(c)$ & $\rhotilde$ & $\rhopw$ & $\wbar$ (finite rule) & $\E[1/\Jinfo]$ (finite rule) \\", r"\midrule"]
for r in gl:
    c = float(r["c"]); A = float(r["mean_info"]); ms = r["msem"]
    msem = r"$\infty$" if ms in ("Inf", "inf") else (sci(float(ms)) if float(ms) > 1e4 else fmt(float(ms), 3))
    wb = float(r["w_bar_truncated"]); wbs = fmt(wb, 4) if wb > 1e-4 else (sci(wb) if wb > 0 else "0")
    lines.append(f"{c:g} & {A:.1f} & {fmt(float(r['rho_tilde']), 4)} & {fmt(float(r['rho_psd']), 4)} & {wbs} & {msem} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", "", r"\smallskip", r"\begin{minipage}{0.95\textwidth}", r"\footnotesize",
          r"\emph{Note.} Evaluated on 60,011 density-weighted nodes with $v=1$. As the multiplier grows the mean-information index approaches one (\cref{prop:overlap}), the pointwise index falls after its maximum (\cref{prop:vanish}), and the finite-rule inverse-information index collapses because the mean reciprocal information explodes; at $c=500$ the finite rule returns an infinite mean reciprocal. The multipliers beyond 10 lie outside any interval used in this paper.",
          r"\end{minipage}", r"\end{table}"]
(TAB / "tab_guttman_limit.tex").write_text("\n".join(lines) + "\n")

# ---------------------------------------------------------------- tail truncation
tt = read_csv(EXTRA["tail_truncation"]); assert len(tt) == 6
lines = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{The truncation artifact: indices of one form under a standardized $t_5$ and a normal population as the quadrature range widens}",
         r"\label{tab:tail-trunc}", r"\small", r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{rrrrrrrr}", r"\toprule",
         r" & & \multicolumn{3}{c}{standardized $t_5$} & \multicolumn{3}{c}{standard normal} \\", r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
         r"Range $\pm$ & excluded mass & $\rhotilde$ & $\rhopw$ & $\wbar$ (finite rule) & $\rhotilde$ & $\rhopw$ & $\wbar$ \\", r"\midrule"]
for r in tt:
    wb = float(r["t5_w_bar"]); wbs = fmt(wb, 4) if wb > 1e-4 else sci(wb)
    lines.append(f"{float(r['range']):g} & {sci(float(r['excluded_mass']))} & {fmt(float(r['t5_rho_tilde']), 4)} & {fmt(float(r['t5_rho_psd']), 4)} & {wbs} & {fmt(float(r['normal_rho_tilde']), 4)} & {fmt(float(r['normal_rho_psd']), 4) if 'normal_rho_psd' in r else '---'} & {fmt(float(r['normal_w_bar']), 4)} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", "", r"\smallskip", r"\begin{minipage}{0.95\textwidth}", r"\footnotesize",
          r"\emph{Note.} The same 20-item 2PL form is evaluated on density-weighted grids truncated at $\pm$ the stated range, with the represented distribution restandardized to variance one. Under the normal population the three indices are stable in every decimal shown. Under the $t_5$ population the mean-information and pointwise indices are stable while the finite-rule inverse-information index falls toward zero as more of the tail is represented, which is the finite form of \cref{prop:tail}: its population value is zero at every range.",
          r"\end{minipage}", r"\end{table}"]
(TAB / "tab_tail_truncation.tex").write_text("\n".join(lines) + "\n")

# ---------------------------------------------------------------- runtime
rt = read_csv(EXTRA["T_s2_runtime"])
algs = ["EQC", "SAC (rho-tilde)", "SAC (w-bar)", "SAC (superpopulation)"]
lab = {"EQC": r"EQC (0.3.0, fixed form)", "SAC (rho-tilde)": r"SAC, $\rhotilde$ (0.3.0, fixed form)",
       "SAC (w-bar)": r"SAC, $\wbar$ (0.3.0, fixed form)", "SAC (superpopulation)": r"SAC, form population (0.3.1)"}
lines = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{Median elapsed seconds per calibration by algorithm and test length}",
         r"\label{tab:runtime}", r"\small", r"\begin{tabular}{lrrrr}", r"\toprule",
         r"Algorithm & 10 items & 20 items & 40 items & 60 items \\", r"\midrule"]
rt_lookup = {(r["algorithm"], r["n_items"]): r for r in rt}
for a in algs:
    cells = []
    for I in ("10", "20", "40", "60"):
        r = rt_lookup.get((a, I)); cells.append("---" if r is None else f"{float(r['median_sec']):.1f}")
    lines.append(f"{lab[a]} & " + " & ".join(cells) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", "", r"\smallskip", r"\begin{minipage}{0.95\textwidth}", r"\footnotesize",
          r"\emph{Note.} Medians over each algorithm's successful runs on the recorded computing setup (30 workers); between 63 and 191 runs per cell. EQC used 20,000 nodes; the fixed-form stochastic calibrations used a 20,000-node presample, 1,000 iterations and 1,000 fresh draws per iteration; the form-population arm was run only at 20 and 40 items. The fixed iteration budget does not make the stochastic cost independent of length. These are measurements of one implementation, not rankings.",
          r"\end{minipage}", r"\end{table}"]
(TAB / "tab_runtime.tex").write_text("\n".join(lines) + "\n")

# ---------------------------------------------------------------- holdout by cell
ms = [r for r in read_csv(SOURCES["s3_m_sensitivity"]) if r["status"] == "ok"]; assert len(ms) == 1680
cells = defaultdict(list)
for r in ms: cells[(r["model"], r["latent_shape"], r["n_items"], int(r["M"]))].append(float(r["holdout_error"]))
Ms = sorted({int(r["M"]) for r in ms}); model_lab = {"rasch": "Rasch", "2pl": "2PL", "3pl_g20": "3PL, $g=.20$"}
lines = [r"\begin{table}[!htbp]", r"\centering",
         r"\caption{Fresh-draw mean absolute error of the 12 holdout cells at each node count (20 seeds per cell)}",
         r"\label{tab:holdout-cell}", r"\footnotesize", r"\setlength{\tabcolsep}{2.2pt}",
         r"\begin{tabular}{lll" + "r" * len(Ms) + "}", r"\toprule",
         r"Model & Shape & $I$ & " + " & ".join(("$M$ = " if M == 500 else "") + (f"{M/1000:g}k" if M >= 1000 else str(M)) for M in Ms) + r" \\", r"\midrule"]
rows_csv = []
for m in ("rasch", "2pl", "3pl_g20"):
    for sh in ("normal", "skew_pos"):
        for I in ("10", "40"):
            vals = [st.mean(abs(x) for x in cells[(m, sh, I, M)]) for M in Ms]
            assert all(len(cells[(m, sh, I, M)]) == 20 for M in Ms)
            lines.append(f"{model_lab[m]} & {shape_lab[sh]} & {I} & " + " & ".join(fmt(v, 4) for v in vals) + r" \\")
            rows_csv.append({"model": m, "shape": sh, "n_items": I, **{f"M{M}": v for M, v in zip(Ms, vals)}})
lines += [r"\bottomrule", r"\end{tabular}", "", r"\smallskip", r"\begin{minipage}{0.95\textwidth}", r"\footnotesize",
          r"\emph{Note.} Target .70; each entry is the mean absolute error, over 20 seeds, of the calibrated form's mean-information index on 200,000 fresh draws (their sample variance) against the target. The pooled means over the 12 cells are the fresh-draw line of \cref{fig:holdout}. Within a node count the cells differ because information is averaged over different populations and forms; the skewed population is the harder case at every node count, whereas test length has no consistent ordering across node counts.",
          r"\end{minipage}", r"\end{table}"]
(TAB / "tab_holdout_by_cell.tex").write_text("\n".join(lines) + "\n")
write_csv(DERIVED / "tab_holdout_by_cell.csv", rows_csv)

# ---------------------------------------------------------------- Arm A, all conditions
T = read_csv(SOURCES["T_s6_arm_a"])
rows = {(r["method"], r["target_rho"], r["te_mean"], r["te_sd"]): r for r in T if r["te_scale"] == "theta_fixed"}
targets = ["0.5", "0.6", "0.7", "0.8", "0.887", "0.9"]
conds = [("0", "0"), ("0", "0.2"), ("0", "0.4"), ("0.4", "0"), ("0.4", "0.2"), ("0.4", "0.4")]
lines = [r"\begingroup\footnotesize\setstretch{1.0}\setlength{\tabcolsep}{3.5pt}", r"\begin{longtable}{llrrrr}",
         r"\caption{All 36 conditions of the fixed-length illustration under the ability-scale convention: rejection rates of the constant- and varying-effect fits (1,000 replications each)}\label{tab:arm-a-all}\\",
         r"\toprule", r"Reference mean & item-effect SD & target & constant effect & varying effect & 95\% interval (constant) \\", r"\midrule", r"\endfirsthead",
         r"\toprule", r"Reference mean & item-effect SD & target & constant effect & varying effect & 95\% interval (constant) \\", r"\midrule", r"\endhead",
         r"\midrule\multicolumn{6}{r}{\footnotesize continued on the next page}\\", r"\endfoot", r"\bottomrule", r"\endlastfoot"]
rows_csv = []
for (m, s) in conds:
    for t in targets:
        rc = rows[("Constant", t, m, s)]; rv = rows[("Varying", t, m, s)]
        n = int(rc["n_reps"]); assert n == 1000 and int(rv["n_reps"]) == 1000
        pc = float(rc["reject_rate"]); pv = float(rv["reject_rate"]); x = round(pc * n)
        lo, hi = clopper_pearson(x, n)
        ci = f"[{fmt(lo)}, {fmt(hi)}]" if (m, s) in (("0.4", "0"), ("0", "0.4")) else ""
        lines.append(f"{fmt(float(m),1) if float(m) else '0'} & {fmt(float(s),1) if float(s) else '0'} & {fmt(float(t), 3).rstrip('0')} & {fmt(pc)} & {fmt(pv)} & {ci} \\\\")
        rows_csv.append({"te_mean": m, "te_sd": s, "target_rho": t, "n_reps": n, "reject_constant": pc, "reject_varying": pv,
                         "cp_lower_constant": lo, "cp_upper_constant": hi, "se_ratio_constant": rc["se_ratio"], "se_ratio_varying": rv["se_ratio"]})
lines += [r"\end{longtable}", r"\endgroup"]
(TAB / "tab_arm_a_all.tex").write_text("\n".join(lines) + "\n")
write_csv(DERIVED / "tab_arm_a_all.csv", rows_csv)

provenance("tables_v5", {**EXTRA, "s3_m_sensitivity": SOURCES["s3_m_sensitivity"], "T_s6_arm_a": SOURCES["T_s6_arm_a"]},
           [TAB / n for n in ("tab_turning_points.tex", "tab_guttman_limit.tex", "tab_tail_truncation.tex", "tab_runtime.tex",
                               "tab_holdout_by_cell.tex", "tab_arm_a_all.tex")],
           {"turning_points": {"psd_min": min(tp_psd), "psd_max": max(tp_psd), "wbar_min": min(tp_w), "wbar_max": max(tp_w)}})
print("turning points: psd", round(min(tp_psd), 2), round(max(tp_psd), 2), "wbar", round(min(tp_w), 2), round(max(tp_w), 2))
print("wrote 6 tables")
