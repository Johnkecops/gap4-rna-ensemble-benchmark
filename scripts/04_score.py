"""Stage E (Gap 4): p(unpaired) vs reactivity, base-pair F1/MCC, runtime/RAM vs length.
Outputs results/scores_per_rna.tsv, results/leaderboard.tsv, results/scaling.tsv, results/fig_scaling.png"""
import os, glob, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C = json.load(open("config.json")); rng = np.random.default_rng(C["split_seed"])
S = pd.read_csv("data/processed/split.tsv", sep="\t").set_index("rna_id")
R = pd.read_csv("data/processed/rnas.tsv", sep="\t").set_index("rna_id")
X = pd.read_csv("data/processed/reactivity.tsv.gz", sep="\t")


def pairs_db(db):
    st, P = [], set()
    for i, c in enumerate(db, 1):
        if c == "(": st.append(i)
        elif c == ")" and st: P.add((st.pop(), i))
    return P


rows, scal = [], []
for f in sorted(glob.glob("results/ensembles/*/*.npz")):
    rid = f.split("/")[-2]; tool = os.path.basename(f)[:-4]
    if rid not in S.index: continue
    z = np.load(f); pu = z["punp"]
    scal.append(dict(rna_id=rid, tool=tool, length=int(z["length"]), runtime_s=float(z["runtime_s"]), peak_rss_mb=float(z["peak_rss_mb"])))
    if S.loc[rid, "set"] != "held_out": continue          # score on held-out RNAs only
    prof = X[X.rna_id == rid]
    if S.loc[rid, "source"] == "MERGE-RNA":               # Rep1 is MERGE fit data -> score on Rep2
        prof = prof[prof.profile.str.endswith("|2")]
        mf = f"results/ensembles/{rid}/merge_shared_b-1.npz"   # same trimmed window for every tool
        if os.path.exists(mf):
            win = np.flatnonzero(~np.isnan(np.load(mf)["punp"])) + 1; prof = prof[prof.pos.isin(win)]
    sp, au = [], []
    for _, g in prof.groupby("profile"):
        p = pu[g.pos.values - 1]; r = g.reactivity.values
        if len(r) < 10 or np.std(p) == 0: continue
        sp.append(spearmanr(p, r).statistic)
        y = r > np.quantile(r, C["reactive_quantile"])
        if 0 < y.sum() < len(y): au.append(roc_auc_score(y, p))
    ref = R.loc[rid, "ref_structure"] if isinstance(R.loc[rid, "ref_structure"], str) else ""
    f1 = mcc = np.nan
    if ref:
        T = pairs_db(ref); Pp = {tuple(x) for x, q in zip(z["ij"].tolist(), z["p"]) if q > 0.5}
        tp = len(T & Pp); ppv = tp / len(Pp) if Pp else 0; sen = tp / len(T) if T else 0
        f1 = 2 * ppv * sen / (ppv + sen) if ppv + sen else 0; mcc = np.sqrt(ppv * sen)
    rows.append(dict(rna_id=rid, tool=tool, length=len(pu), stratum=S.loc[rid, "stratum"], source=S.loc[rid, "source"],
                     spearman=np.nanmean(sp) if sp else np.nan, auc=np.nanmean(au) if au else np.nan, bp_f1=f1, bp_mcc=mcc))

D = pd.DataFrame(rows); D.to_csv("results/scores_per_rna.tsv", sep="\t", index=False)
SC = pd.DataFrame(scal); SC.to_csv("results/runtime_memory.tsv", sep="\t", index=False)


def boot(v):
    v = v[~np.isnan(v)]
    if len(v) < 2: return (np.nan, np.nan)
    b = [rng.choice(v, len(v)).mean() for _ in range(C["bootstrap"])]
    return tuple(np.percentile(b, [2.5, 97.5]))


L = []
for (st, tool), g in D.groupby(["stratum", "tool"]):
    lo, hi = boot(g.spearman.values); alo, ahi = boot(g.auc.values)
    L.append(dict(stratum=st, tool=tool, n_rna=g.spearman.notna().sum(), spearman=g.spearman.mean(), sp_lo=lo, sp_hi=hi,
                  auc=g.auc.mean(), auc_lo=alo, auc_hi=ahi, n_ref=g.bp_f1.notna().sum(), bp_f1=g.bp_f1.mean(), bp_mcc=g.bp_mcc.mean()))
LB = pd.DataFrame(L).sort_values(["stratum", "spearman"], ascending=[True, False])
LB.to_csv("results/leaderboard.tsv", sep="\t", index=False, float_format="%.3f")

# scaling: log-log slope of runtime / memory vs length
sl = []
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for tool, g in SC.groupby("tool"):
    g = g[g.length >= 50]
    if len(g) < 3: continue
    a = np.polyfit(np.log10(g.length), np.log10(g.runtime_s.clip(1e-4)), 1)[0]
    b = np.polyfit(np.log10(g.length), np.log10(g.peak_rss_mb.clip(1)), 1)[0]
    sl.append(dict(tool=tool, n=len(g), max_len=g.length.max(), runtime_slope=a, memory_slope=b,
                   median_runtime_s=g.runtime_s.median(), max_rss_mb=g.peak_rss_mb.max()))
    ax[0].loglog(g.length, g.runtime_s, ".", label=tool, ms=3); ax[1].loglog(g.length, g.peak_rss_mb, ".", ms=3)
ax[0].set(xlabel="length (nt)", ylabel="runtime (s)"); ax[1].set(xlabel="length (nt)", ylabel="peak RSS (MB)")
ax[0].legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig("results/fig_scaling.png", dpi=150)
pd.DataFrame(sl).to_csv("results/scaling.tsv", sep="\t", index=False, float_format="%.3f")
print(LB.to_string(index=False)); print(pd.DataFrame(sl).to_string(index=False))
