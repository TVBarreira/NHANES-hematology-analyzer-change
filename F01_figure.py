"""F01_figure.py - draws the Research Letter figure from M02's outputs.

Reads results/m_all_boundary_steps.csv (panel A) and results/m_figure_panels.csv
(panel B). Run from the project folder, or set NHANES_BIOAGE_DIR.
Requires Python 3 with matplotlib and numpy. Writes results/figures/Figure.png.

Panel B uses the corrected series expressed against the observed age fit
("RDW + WBC (same reference)"), so the two series coincide before 2013-2014.
If an older m_figure_panels.csv lacks that series, the constant pre-2013 offset
between the two residual series is removed instead, which is numerically the
same thing.
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

base = os.environ.get("NHANES_BIOAGE_DIR", os.getcwd())
R = os.path.join(base, "results")
os.makedirs(os.path.join(R, "figures"), exist_ok=True)
steps = list(csv.DictReader(open(os.path.join(R, "m_all_boundary_steps.csv"), encoding="utf-8-sig")))
panels = list(csv.DictReader(open(os.path.join(R, "m_figure_panels.csv"), encoding="utf-8-sig")))

ROWS = [("rdw", "Red cell distribution width"), ("wbc", "White blood cell count"), ("rbc", "Red blood cell count"),
        ("hct", "Hematocrit"), ("mch", "Mean corpuscular hemoglobin"), ("mchc", "MCH concentration"),
        ("mono", "Monocytes, %"), ("hgb", "Hemoglobin"), ("mcv", "Mean corpuscular volume"),
        ("lymph", "Lymphocytes, %"), ("neut", "Neutrophils, %"), ("plt", "Platelets")]
ROWS = [(a, n) for a, n in ROWS if any(s["analyte"] == a for s in steps)]
YRS = [1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2022]
LAB = ["1999-\n2000", "2001-\n2002", "2003-\n2004", "2005-\n2006", "2007-\n2008", "2009-\n2010",
       "2011-\n2012", "2013-\n2014", "2015-\n2016", "2017-\n2018", "2021-\n2023"]
RED, BLUE, GREY = "#B00020", "#1F4E79", "#9AA3AD"

fig = plt.figure(figsize=(7.2, 10.6), dpi=300)
gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.5)

# ---- A: every measure, change at the analyzer boundary against the other 9 ----
ax = fig.add_subplot(gs[0, 0]); n = len(ROWS)
for i, (a, name) in enumerate(ROWS):
    y = n - 1 - i
    r = [x for x in steps if x["analyte"] == a]
    other = [float(x["pct"]) for x in r if x["boundary"] != "2011-2013"]
    gh = [x for x in r if x["boundary"] == "2011-2013"][0]
    v, lo, hi, p = (float(gh[k]) for k in ("pct", "pct_lo", "pct_hi", "p"))
    if i % 2 == 0: ax.axhspan(y - .5, y + .5, color="#F5F6F8", zorder=0)
    ax.scatter(other, [y] * len(other), s=16, color=GREY, zorder=2, linewidths=0)
    ax.plot([lo, hi], [y, y], color=RED, lw=1.4, zorder=3)
    ax.scatter([v], [y], s=46, zorder=4, color=RED if p < .05 else "white", edgecolors=RED, linewidths=1.4)
ax.axvline(0, color="#555", lw=.8, zorder=1)
ax.set_yticks(range(n)); ax.set_yticklabels([nm for _, nm in ROWS][::-1], fontsize=9); ax.set_ylim(-.6, n - .4)
ax.set_xlabel("Change between adjacent cycles, % (adolescents aged 12-19 years)", fontsize=9)
ax.tick_params(axis="x", labelsize=8.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.set_title("A  Blood-count measures", fontsize=11, loc="left", fontweight="bold")
ax.legend(handles=[
    Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY, markeredgecolor="none", markersize=5, label="Other 9 cycle boundaries"),
    Line2D([0], [0], marker="o", color=RED, markerfacecolor=RED, markersize=6.5, lw=1.4, label="Analyzer change, 2011-12 to 2013-14 (95% CI)"),
    Line2D([0], [0], marker="o", color=RED, markerfacecolor="white", markersize=6.5, lw=1.4, label="Analyzer change, not significant")],
    fontsize=7.8, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=2)

# ---- B: phenotypic age acceleration by cycle, measured vs corrected ----
def series(tag):
    d = {int(r["year"]): r for r in panels if r["panel"] == "C" and r["series"] == tag}
    if len(d) < len(YRS): return None
    return tuple(np.array([float(d[y][k]) for y in YRS]) for k in ("mean", "lo", "hi"))
mo, lo_o, hi_o = series("none")
same = series("RDW + WBC (same reference)")
if same is None:                                   # older M02 output: remove the constant pre-2013 offset
    mc, lo_c, hi_c = series("RDW + WBC"); off = float(np.mean(mc[:7] - mo[:7]))
    mc, lo_c, hi_c = mc - off, lo_c - off, hi_c - off
else:
    mc, lo_c, hi_c = same
x = np.arange(len(YRS))
bx = fig.add_subplot(gs[1, 0])
bx.axvspan(6.5, 10.5, color="#FBE9EC", zorder=0); bx.axvline(6.5, color="#444", lw=1.2, ls="--", zorder=1)
bx.fill_between(x, lo_o, hi_o, color=RED, alpha=.13, zorder=2); bx.plot(x, mo, marker="o", color=RED, lw=1.8, ms=4.5, label="As measured", zorder=4)
bx.fill_between(x[6:], lo_c[6:], hi_c[6:], color=BLUE, alpha=.15, zorder=2)
bx.plot(x[6:], mc[6:], marker="s", color=BLUE, lw=1.8, ms=4.5, label="After correcting RDW and WBC", zorder=3)
bx.plot([6, 7.45], [mo[6], mo[6]], color=RED, lw=.8, ls=":", zorder=3)
bx.annotate("", xy=(7.45, mo[7]), xytext=(7.45, mo[6]), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.3))
bx.text(7.6, (mo[6] + mo[7]) / 2 + .25, f"{mo[7]-mo[6]:+.1f} years", color=RED, fontsize=9.5, fontweight="bold", va="center")
bx.text(7.6, mc[7] - .35, f"{mc[7]-mc[6]:+.1f} years", color=BLUE, fontsize=9, fontweight="bold", va="center")
yl = bx.get_ylim(); bx.text(6.6, yl[1] - .05 * (yl[1] - yl[0]), "New analyzer", fontsize=8, color="#444", va="top")
bx.set_xticks(x); bx.set_xticklabels(LAB, fontsize=7.6)
bx.set_ylabel("Mean phenotypic age acceleration, years", fontsize=9); bx.tick_params(axis="y", labelsize=8.5)
bx.legend(fontsize=8, frameon=False, loc="upper left")
for sp in ("top", "right"): bx.spines[sp].set_visible(False)
bx.set_title("B  Phenotypic age acceleration (adults)", fontsize=11, loc="left", fontweight="bold")
out = os.path.join(R, "figures", "Figure.png")
plt.savefig(out, bbox_inches="tight", facecolor="white"); print("saved", out)
