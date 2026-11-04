"""Positive control (doc: Datasets): ViennaRNA 2.7.2 on EternaBench Cloud Lab Round 1 demo set,
scored with EternaBench's own concatenation + 95th-pct filter + polyA exclusion, Pearson.
Published vienna_2 value: data/EternaBench_ChemMapping_Example_Correlations_Round1.csv."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import pearsonr
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EB = f"{ROOT}/data/raw/EternaBench"; os.environ.setdefault("ETERNABENCH_PATH", EB); sys.path.insert(0, EB); sys.path.insert(0, f"{ROOT}/scripts")
import foldtools as ft

d = pd.read_json(f"{EB}/data/DEMO_ChemMapping.json.zip")
d = d[d.passed_CDHIT_filter].reset_index(drop=True)
d["p_vienna_2"] = [list(ft.vienna(s, k=0)[1][[x for x in sp if x < len(s)]]) for s, sp in zip(d.sequence, d.seqpos)]
# same steps as eternabench.chemmapping_utils (write_concatenated_dataframe, get_polyA_indicator(6), filter_data(95));
# re-implemented because that module imports the retired `RDATKit.rdatkit` package name
rows = []
for s, sp, r, p in zip(d.sequence, d.seqpos, d.reactivity, d.p_vienna_2):
    t = [s[x] for x in sp if x < len(s)]
    pa = [0] * 5 + [int("".join(t[i - 6:i]) == "A" * 6) for i in range(6, len(t) + 1)]
    rows += zip(r[:len(t)], p, pa)
c = pd.DataFrame(rows, columns=["reactivity", "p_vienna_2", "in_polyA"]).dropna()
c = c[c.reactivity > 0]; c = c[c.reactivity < np.percentile(c.reactivity, 95)]; c = c[c.in_polyA == 0]
r = pearsonr(c.p_vienna_2, c.reactivity).statistic
pub = pd.read_csv(f"{EB}/data/EternaBench_ChemMapping_Example_Correlations_Round1.csv").set_index("package").loc["vienna_2", "mean"]
out = f"constructs={len(d)} nucleotides={len(c)} pearson_ours={r:.4f} published_vienna_2={pub:.4f} diff={r - pub:+.4f}"
print(out); open(f"{ROOT}/results/positive_control_eternabench.txt", "w").write(out + "\n")
