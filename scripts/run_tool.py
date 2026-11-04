"""Stage D worker: python run_tool.py RNA_ID TOOL BEAM -> results/ensembles/RNA_ID/TOOL_bBEAM.npz
TOOL in {vienna, linearpartition, linearsampling}. Stores bpp (i,j,p), punp, samples,
runtime_s, peak_rss_mb (GNU time for binaries; RSS increment over baseline for in-process ViennaRNA)."""
import sys, os, time, json, resource, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foldtools as ft
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rid, tool, beam = sys.argv[1], sys.argv[2], int(sys.argv[3])
C = json.load(open("config.json"))
R = pd.read_csv("data/processed/rnas.tsv", sep="\t").set_index("rna_id")
seq = R.loc[rid, "seq"]; n = len(seq)
out = f"results/ensembles/{rid}/{tool}_b{beam}.npz"; os.makedirs(os.path.dirname(out), exist_ok=True)
if (tool == "vienna" or beam == 0) and n > C["exact_max_len"]:
    print("skip: exact run above exact_max_len"); sys.exit(0)
base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
t0 = time.time(); samples = []
if tool == "vienna":
    bpp, pu, samples = ft.vienna(seq, C["n_samples"], C["bpp_cutoff"]); who = "self"
elif tool == "linearpartition":
    bpp, pu = ft.linearpartition(seq, beam, C["bpp_cutoff"]); who = "child"
elif tool == "linearsampling":
    samples = ft.linearsampling(seq, beam, C["n_samples"]); who = "child"
    pu = ft.punp_from_samples(samples, n)
    bpp = {}
    for s in samples:   # sample pair frequencies
        st = []
        for i, c in enumerate(s, 1):
            if c == "(": st.append(i)
            elif c == ")": j = st.pop(); bpp[(j, i)] = bpp.get((j, i), 0) + 1 / len(samples)
dt = time.time() - t0
rss = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 - base) if who == "self" else ft.LAST_RSS_MB[0]
ij = np.array(list(bpp.keys()), dtype=np.int32).reshape(-1, 2); p = np.array(list(bpp.values()), dtype=np.float32)
np.savez_compressed(out, ij=ij, p=p, punp=pu.astype(np.float32), samples=np.array(samples),
                    runtime_s=dt, peak_rss_mb=rss, length=n)
print(f"{rid} {tool} b{beam} n={n} {dt:.2f}s {rss:.0f}MB")
