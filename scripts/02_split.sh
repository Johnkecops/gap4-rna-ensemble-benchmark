#!/bin/bash
# Stage C: cd-hit-est clustering -> whole clusters to MERGE-RNA fit set or held-out set;
# length strata. Output: data/processed/split.tsv (rna_id length stratum cluster set)
set -euo pipefail
cd "$(dirname "$0")/.."
P=data/processed
ID=$(python -c "import json;print(json.load(open('config.json'))['cdhit_identity'])")
W=$(python -c "import json;print(json.load(open('config.json'))['cdhit_word'])")
python - <<'EOF'
import pandas as pd
R = pd.read_csv("data/processed/rnas.tsv", sep="\t")
with open("data/processed/all.fa", "w") as f:
    for r in R.itertuples(): f.write(f">{r.rna_id}\n{r.seq}\n")
EOF
cd-hit-est -i $P/all.fa -o $P/clusters -c $ID -n $W -M 4000 -T 4 -d 0 > logs/cdhit.log
python - <<'EOF'
import json, random, pandas as pd
C = json.load(open("config.json")); random.seed(C["split_seed"])
cl, cur = {}, None
for line in open("data/processed/clusters.clstr"):
    if line.startswith(">Cluster"): cur = int(line.split()[1])
    else: cl[line.split(">")[1].split("...")[0]] = cur
R = pd.read_csv("data/processed/rnas.tsv", sep="\t"); R["cluster"] = R.rna_id.map(cl)
# MERGE-RNA can only be fitted on MaP-count data (source MERGE-RNA); split those clusters
mclus = sorted(R.loc[R.source == "MERGE-RNA", "cluster"].unique()); random.shuffle(mclus)
fit = set(mclus[: max(1, round(len(mclus) * C["fit_fraction"]))])
R["set"] = ["merge_fit" if c in fit else "held_out" for c in R.cluster]
s = C["length_strata"]
R["stratum"] = pd.cut(R.length, [0, s[0], s[1], s[2], 10**9], right=False,
                      labels=[f"<{s[0]}", f"{s[0]}-{s[1]}", f"{s[1]}-{s[2]}", f">{s[2]}"])
R[["rna_id", "length", "stratum", "cluster", "set", "source"]].to_csv("data/processed/split.tsv", sep="\t", index=False)
print(R.groupby(["set", "stratum"], observed=True).size())
EOF
