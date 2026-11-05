"""Stage G (figures): secondary-structure visualisation for samples fixed in config.json["viz"].
 fig_qc_misaligned.png        profiles rejected because positions and values could not be aligned (QC diagnostic)
 fig_reactivity_on_structure.png  normalised reactivity on the reference structure (aligned profiles)
 fig_punp_vs_reactivity.png   exact p(unpaired), beam-100 p(unpaired) and reactivity on held-out RNAs
Reads existing files only; nothing is refolded."""
import os, json, numpy as np, pandas as pd, RNA
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from rdat_kit import RDATFile
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C = json.load(open("config.json"))["viz"]
R = pd.read_csv("data/processed/rnas.tsv", sep="\t").set_index("rna_id")
X = pd.read_csv("data/processed/reactivity.tsv.gz", sep="\t")
SC = pd.read_csv("results/scores_per_rna.tsv", sep="\t").set_index(["rna_id", "tool"])
CMAP = plt.get_cmap("viridis")


def pairs(db):
    st, P = [], []
    for i, c in enumerate(db):
        if c == "(": st.append(i)
        elif c == ")": P.append((st.pop(), i))
    return P


def draw_structure(ax, db, values, title, vmax=1.0, label=""):
    """naview layout; nodes coloured by values (NaN = grey)."""
    co = RNA.naview_xy_coordinates(db); n = len(db)
    x = np.array([co[i].X for i in range(n)]); y = np.array([co[i].Y for i in range(n)])
    ax.plot(x, y, "-", color="0.6", lw=0.6, zorder=1)
    for i, j in pairs(db): ax.plot([x[i], x[j]], [y[i], y[j]], "-", color="0.35", lw=0.6, zorder=1)
    v = np.asarray(values, float)
    ok = ~np.isnan(v)
    ax.scatter(x[~ok], y[~ok], s=9, c="0.85", edgecolors="0.6", linewidths=0.3, zorder=2)
    sc = ax.scatter(x[ok], y[ok], s=9, c=np.clip(v[ok], 0, vmax), cmap=CMAP, vmin=0, vmax=vmax, edgecolors="none", zorder=3)
    ax.text(x[0], y[0], " 5'", fontsize=6, ha="left", va="center")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=7.5)
    cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.01); cb.ax.tick_params(labelsize=6); cb.set_label(label, fontsize=6.5)


def react(rid, rep=1):
    g = X[X.profile == X[X.rna_id == rid].profile.unique()[rep - 1]]
    v = np.full(R.loc[rid, "length"], np.nan); v[g.pos.values - 1] = g.reactivity.values
    return v, g.profile.iloc[0].split("|")[1]


def ens(rid, tool):
    return np.load(f"results/ensembles/{rid}/{tool}.npz")["punp"].astype(float)


# ---------------- Fig. A: QC diagnostic for unalignable profiles
fig, axes = plt.subplots(len(C["qc_mismatch"]), 1, figsize=(7.2, 2.3 * len(C["qc_mismatch"])))
for ax, rid in zip(np.atleast_1d(axes), C["qc_mismatch"]):
    r = RDATFile(); r.load(open(f"data/raw/rmdb/{rid}.rdat"))
    def qc_lines(c):   # same selection as 01_qc.py: WT, no ligand, first two, reactivity datatypes, not Nomod
        wt = [d for d in c.data if (d.annotations.get("mutation") or ["WT"])[0].upper() in ("WT", "") and not d.annotations.get("ligand")][:2]
        return [d for d in wt if (d.annotations.get("modifier") or ["x"])[0].lower() not in ("nomod", "none", "?")
                and (d.annotations.get("datatype") or ["REACTIVITY"])[0].upper().startswith("REACTIVITY")]
    c, d = next((c, d) for c in r.constructs.values() for d in qc_lines(c)   # the line QC rejected
                if len(d.seqpos or c.seqpos) != len(d.values))
    vals = np.array(d.values, float); sp = np.array(d.seqpos if d.seqpos else c.seqpos)
    k = np.arange(1, len(vals) + 1)
    ax.axvspan(0.5, len(sp) + 0.5, color="#4a78b5", alpha=0.08, label=f"SEQPOS positions of the full construct ({len(sp)})")
    ax.bar(k, vals, width=1.0, color="#d9534f", label=f"values in the data line ({len(vals)})")
    ax.set_xlim(0, max(len(vals), len(sp)) + 1); ax.set_xlabel("index", fontsize=7.5)
    ax.set_ylabel("raw reactivity", fontsize=7.5); ax.tick_params(labelsize=6.5)
    own = (d.annotations.get("sequence") or [""])[0]
    why = (f"line carries its own {len(own)}-nt sequence ({d.annotations.get('ID', [''])[0]}): a transcript-length fragment, "
           "not the full-length RNA") if own else "fewer values than positions: mapping undefined"
    ax.set_title(f"{rid}: construct {len(c.sequence)} nt, {len(sp)} SEQPOS positions, {len(vals)} values, "
                 f"{(d.annotations.get('modifier') or [''])[0]}\n{why}; excluded", fontsize=7.5)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig("results/fig_qc_misaligned.png", dpi=300); plt.close(fig)

# ---------------- Fig. B: reactivity on the reference structure
fig, axes = plt.subplots(1, len(C["aligned"]), figsize=(7.2, 3.8))
for ax, rid in zip(axes, C["aligned"]):
    v, probe = react(rid, 1)
    draw_structure(ax, R.loc[rid, "ref_structure"], v, f"{rid} ({R.loc[rid, 'length']} nt, {probe})\nreference structure",
                   vmax=2.0, label="normalised reactivity")
fig.tight_layout(); fig.savefig("results/fig_reactivity_on_structure.png", dpi=300); plt.close(fig)

# ---------------- Fig. C: p(unpaired) vs reactivity on held-out RNAs
rows = C["heldout"]; lin = C["linear_tool"]
fig = plt.figure(figsize=(7.2, 2.6 * len(rows)))
gs = fig.add_gridspec(len(rows), 3, width_ratios=[1, 1, 1.6], hspace=0.55, wspace=0.25)
for r_i, rid in enumerate(rows):
    db = R.loc[rid, "ref_structure"] if isinstance(R.loc[rid, "ref_structure"], str) else ""; n = R.loc[rid, "length"]
    pe = ens(rid, "vienna_b-1"); pl = ens(rid, lin); v, probe = react(rid, 1)
    rho_e = SC.loc[(rid, "vienna_b-1"), "spearman"]; rho_l = SC.loc[(rid, lin), "spearman"]
    head = f"{rid} ({n} nt, {probe})"
    if db and n <= C["max_2d_len"]:
        draw_structure(fig.add_subplot(gs[r_i, 0]), db, pe, head + "\nexact p(unpaired)", label="p(unpaired)")
        draw_structure(fig.add_subplot(gs[r_i, 1]), db, v, "reactivity", vmax=2.0, label="norm. reactivity")
    else:   # arc diagram: exact pairs p>0.5 above, linear-time pairs p>0.5 below
        ax = fig.add_subplot(gs[r_i, 0:2])
        for tool, (t1, t2), col in (("vienna_b-1", (0, 180), "#1f4e8c"), (lin, (180, 360), "#d9534f")):
            z = np.load(f"results/ensembles/{rid}/{tool}.npz")
            for (i, j), p in zip(z["ij"], z["p"]):
                if p > 0.5: ax.add_patch(Arc(((i + j) / 2 - 1, 0), j - i, (j - i) * 0.5, theta1=t1, theta2=t2, lw=0.3, color=col))
        ax.set_xlim(0, n); ax.set_ylim(-n * 0.27, n * 0.27); ax.axhline(0, color="0.5", lw=0.4)
        ax.set_yticks([]); ax.tick_params(labelsize=6); ax.set_xlabel("position", fontsize=7)
        ax.set_title(head + f"\npairs with p > 0.5: exact (above), {lin.replace('_b', ' b')} (below)", fontsize=7.5)
    ax = fig.add_subplot(gs[r_i, 2]); pos = np.arange(1, n + 1); m = ~np.isnan(v)
    ax.bar(pos[m], np.clip(v[m], 0, 2) / 2, width=1.0, color="0.75", label="reactivity / 2 (clipped)")
    ax.plot(pos, pe, color="#1f4e8c", lw=0.7, label="exact p(unpaired)")
    ax.plot(pos, pl, color="#d9534f", lw=0.7, ls="--", label=f"{lin.replace('_b', ' b')} p(unpaired)")
    ax.set_xlim(1, n); ax.set_ylim(0, 1.05); ax.tick_params(labelsize=6); ax.set_xlabel("position", fontsize=7)
    mad = np.mean(np.abs(pe - pl))
    ax.set_title(f"ρ exact {rho_e:.2f}, ρ {lin.split('_')[1]} {rho_l:.2f}, MAD {mad:.3f}", fontsize=7.5)
    if r_i == 0: ax.legend(fontsize=6, frameon=False, loc="upper right")
fig.savefig("results/fig_punp_vs_reactivity.png", dpi=300, bbox_inches="tight"); plt.close(fig)
print("figures written")
