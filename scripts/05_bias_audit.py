"""Stage F (Gap 3a): exact ViennaRNA vs LinearPartition-V / LinearSampling beams on RNAs where
the exact run exists. Outputs results/bias_per_rna.tsv, results/bias_summary.tsv"""
import os, glob, json, numpy as np, pandas as pd
from collections import Counter
from scipy.stats import kendalltau
from scipy.spatial.distance import jensenshannon
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
C = json.load(open("config.json"))
S = pd.read_csv("data/processed/split.tsv", sep="\t").set_index("rna_id")
SC = pd.read_csv("results/scores_per_rna.tsv", sep="\t")


def pdict(z): return {tuple(k): float(v) for k, v in zip(z["ij"].tolist(), z["p"])}


rows = []
for vf in sorted(glob.glob("results/ensembles/*/vienna_b-1.npz")):
    rid = vf.split("/")[-2]; V = np.load(vf); vp = pdict(V)
    vs = Counter(V["samples"].tolist())
    for f in glob.glob(f"results/ensembles/{rid}/linear*_b*.npz"):
        tool = os.path.basename(f)[:-4]; Z = np.load(f); zp = pdict(Z)
        keys = set(vp) | set(zp)
        mad_pair = np.mean([abs(vp.get(k, 0) - zp.get(k, 0)) for k in keys]) if keys else 0
        mad_unp = float(np.mean(np.abs(V["punp"] - Z["punp"])))
        jsd = np.nan
        if Z["samples"].size and int(V["length"]) <= C["jsd_max_len"]:
            zs = Counter(Z["samples"].tolist()); U = sorted(set(vs) | set(zs))
            jsd = jensenshannon([vs.get(u, 0) for u in U], [zs.get(u, 0) for u in U], base=2) ** 2
        rows.append(dict(rna_id=rid, length=int(V["length"]), stratum=S.loc[rid, "stratum"], tool=tool,
                         mad_pair=mad_pair, mad_unpaired=mad_unp, max_abs_unpaired=float(np.max(np.abs(V["punp"] - Z["punp"]))), jsd=jsd))
B = pd.DataFrame(rows)
# change in reactivity correlation + Kendall tau of per-RNA ranking (exact vs linear-time)
ex = SC[SC.tool == "vienna_b-1"].set_index("rna_id").spearman
d_sp, tau = [], {}
for tool, g in SC[SC.tool.str.startswith("linear")].groupby("tool"):
    g = g.set_index("rna_id").spearman; common = g.index.intersection(ex.index).dropna()
    a, b = ex[common].values, g[common].values; ok = ~(np.isnan(a) | np.isnan(b))
    tau[tool] = kendalltau(a[ok], b[ok]).statistic if ok.sum() > 2 else np.nan
    for r in common: d_sp.append(dict(rna_id=r, tool=tool, d_spearman=g[r] - ex[r]))
B = B.merge(pd.DataFrame(d_sp), on=["rna_id", "tool"], how="left") if d_sp else B
B.to_csv("results/bias_per_rna.tsv", sep="\t", index=False, float_format="%.5f")
agg = B.groupby(["tool", "stratum"]).agg(n=("rna_id", "size"), max_len=("length", "max"),
        mad_unp_mean=("mad_unpaired", "mean"), mad_unp_p95=("mad_unpaired", lambda x: np.percentile(x, 95)),
        mad_pair_mean=("mad_pair", "mean"), jsd_mean=("jsd", "mean"),
        d_spearman_mean=("d_spearman", "mean"), d_spearman_maxabs=("d_spearman", lambda x: np.nanmax(np.abs(x)) if x.notna().any() else np.nan)).reset_index()
agg["kendall_tau_rna_rank"] = agg.tool.map(tau)
agg["within_bound"] = agg.mad_unp_p95 <= C["mad_tolerance"]
agg.to_csv("results/bias_summary.tsv", sep="\t", index=False, float_format="%.4f")
print(agg.to_string(index=False))
