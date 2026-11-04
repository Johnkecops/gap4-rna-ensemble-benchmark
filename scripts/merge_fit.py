"""MERGE-RNA in two modes (only MaP-count RNAs shipped with MERGE-RNA can be fitted).
 mode_shared: physical params fitted on 'merge_fit' RNAs (Rep1+Rep2), applied to held-out
              RNAs at 0 mM with no soft constraints -> results/ensembles/<rid>/merge_shared_b-1.npz
 mode_perRNA: per-RNA sequential fit (physical + lambda_sc) on Rep1 only, scored later on Rep2
              -> results/ensembles/<rid>/merge_rep1fit_b-1.npz
Writes p(unpaired)=1-pairing_probs and pairs >cutoff."""
import os, sys, re, glob, json, time, resource, numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MR = f"{ROOT}/tools/MERGE-RNA"; os.environ["RNA_STRUCT_HOME"] = MR; sys.path.insert(0, MR)
import matplotlib; matplotlib.use("Agg")
import RNA
from merge_rna import Experiment, MultiSystemsFit
os.chdir(ROOT)
C = json.load(open("config.json")); S = pd.read_csv("data/processed/split.tsv", sep="\t")
S = S[S.source == "MERGE-RNA"]
mode = sys.argv[1]


def exps(system, reps):
    out = []
    for rep in reps:
        for d in sorted(glob.glob(f"{MR}/data/structured_rnas/{system}/Rep{rep}_*mM")):
            out.append(Experiment(f"{d}/info.txt"))
    return out


FULL = pd.read_csv("data/processed/rnas.tsv", sep="\t").set_index("rna_id").seq


def save(rid, tag, bpp, dt, seq):
    """MERGE-RNA folds the primer-trimmed sequence; embed into full-length coordinates (NaN outside)."""
    off = FULL[rid].find(seq.replace("T", "U")); assert off >= 0, rid
    m = bpp.shape[0]; n = len(FULL[rid]); P = np.triu(bpp, 1); i, j = np.nonzero(P > C["bpp_cutoff"])
    pu_t = 1 - (bpp + bpp.T).sum(0) if np.allclose(np.tril(bpp, -1), 0) else 1 - bpp.sum(0)
    pv = P[i, j]; pu = np.full(n, np.nan); pu[off:off + m] = pu_t; i, j = i + off, j + off
    os.makedirs(f"results/ensembles/{rid}", exist_ok=True)
    np.savez_compressed(f"results/ensembles/{rid}/{tag}_b-1.npz", ij=np.c_[i + 1, j + 1].astype(np.int32),
                        p=pv.astype(np.float32), punp=np.clip(pu, 0, 1).astype(np.float32), samples=np.array([]),
                        runtime_s=dt, peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, length=n)


def fit(train, val, mode_, out, sc):
    m = MultiSystemsFit(experiments=train, validation_exps=val, output_suffix=mode_, root_dir=out,
                        infer_1D_sc=sc, mask_edges=(C["end_mask_nt"], C["end_mask_nt"]), fit_mode=mode_,
                        do_plots=False, print_to_std_out=False, overwrite=True)
    return m, m.fit()


if mode == "shared":
    fitsys = [r.rna_id.replace("MERGE_", "") for r in S.itertuples() if r.set == "merge_fit"]
    t0 = time.time()
    pf = "results/merge_runs/shared/physical_only/detailed/physical_params.txt"
    if not os.path.exists(pf):
        fit([e for s in fitsys for e in exps(s, (1, 2))], None, "physical_only", "results/merge_runs/shared", False)
    txt = open(pf).read()
    mu_r, p_b = (float(re.search(rf"^{k}\s*=\s*(\S+)", txt, re.M).group(1)) for k in ("mu_r", "p_b"))
    print("shared params mu_r", mu_r, "p_b", p_b)
    json.dump({"fit_systems": fitsys, "mu_r": float(mu_r), "p_b": float(p_b)}, open("results/merge_runs/shared_params.json", "w"))
    for r in S[S.set == "held_out"].itertuples():
        e = exps(r.rna_id.replace("MERGE_", ""), (1,))[0]
        kBT = 0.0019872 * (273.15 + 37); beta = 1 / kBT
        mu = mu_r + kBT * np.log(0.1 / 1000)       # 0 mM, as in MERGE-RNA save_pairing_probabilities
        pen = -np.log((1 + np.exp(beta * (mu - p_b))) / (1 + np.exp(beta * mu))) * kBT
        # apply the pairing penalty via MERGE-RNA's own routine
        from merge_rna.fit import ExperimentFit
        ef = ExperimentFit(e, infer_1D_sc=False, eps_b=None, is_training=False, mask_edges=None)
        bpp = ef.get_bpp_matrix(pen, lambda_sc=None)
        save(r.rna_id, "merge_shared", bpp, time.time() - t0, e.seq)
elif mode == "perRNA":
    if len(sys.argv) < 3: sys.exit("perRNA needs a system name (run one process per RNA)")
    only = sys.argv[2]
    for r in S[S.set == "held_out"].itertuples():
        s = r.rna_id.replace("MERGE_", ""); t0 = time.time(); out = f"results/merge_runs/rep1fit_{s}"
        if only and s != only: continue
        if not glob.glob(f"{out}/**/base_pairing_probabilities.txt", recursive=True):
            fit(exps(s, (1,)), None, "sequential", out, True)
        f = glob.glob(f"{out}/**/base_pairing_probabilities.txt", recursive=True)
        if not f: print("no bpp output for", s); continue
        save(r.rna_id, "merge_rep1fit", np.loadtxt(f[0]), time.time() - t0, exps(s, (1,))[0].seq)
        print(s, "done", round(time.time() - t0), "s", flush=True)
