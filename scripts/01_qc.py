"""Stage B: QC + normalise probing profiles -> data/processed/reactivity.tsv.gz
and data/processed/rnas.tsv (rna_id, seq, length, source, probe, temp, ref_structure).
Sources: RMDB RDAT (WT profiles only) and MERGE-RNA DMS-MaP .rc files (natural RNAs)."""
import sys, os, glob, re, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools", "MERGE-RNA"))
from rdat_kit import RDATFile
from merge_rna.experiment import load_rc

CFG = json.load(open("config.json"))
MIN_SNR, END_MASK = CFG["min_snr"], CFG["end_mask_nt"]
rows, rnas, log = [], [], []


def valid_ref(db, n):
    """reference must be dot-bracket of length n with at least one pair and balanced brackets"""
    if not db or len(db) != n or not set(db) <= set("().") or "(" not in db: return False
    depth = 0
    for ch in db:
        depth += (ch == "(") - (ch == ")")
        if depth < 0: return False
    return depth == 0


def norm(v):
    """2-8%-style normalisation: divide by mean of values between 90th and 98th pct."""
    ok = v[np.isfinite(v)]
    if ok.size < 10: return None
    lo, hi = np.percentile(ok, [90, 98]); top = ok[(ok >= lo) & (ok <= hi)]
    s = top.mean() if top.size and top.mean() > 0 else None
    return v / s if s else None


def add(rid, seq, probe, temp, src, pos, vals, ref=None, rep=None):
    seq = seq.upper().replace("T", "U")
    vals = np.asarray(vals, float); pos = np.asarray(pos, int)
    keep = (pos >= 1) & (pos <= len(seq)) & np.isfinite(vals)
    pos, vals = pos[keep], vals[keep]
    if pos.size < 20: log.append((rid, "too_few_positions")); return
    m = (pos > pos.min() + END_MASK - 1) & (pos < pos.max() - END_MASK + 1)   # mask ends
    pos, vals = pos[m], vals[m]
    nt = np.array([seq[p - 1] for p in pos])
    if probe.upper() == "DMS":
        k = np.isin(nt, ["A", "C"]); pos, vals, nt = pos[k], vals[k], nt[k]
    v = norm(vals)
    if v is None: log.append((rid, "normalisation_failed")); return
    prof = f"{rid}|{probe}|{rep or 1}"
    rows.extend(zip([rid] * len(pos), [prof] * len(pos), pos, nt, [probe] * len(pos), v))
    if rid not in {r["rna_id"] for r in rnas}:
        rnas.append(dict(rna_id=rid, seq=seq, length=len(seq), source=src, probe=probe, temp=temp, ref_structure=ref or ""))


# ---- RMDB
for f in sorted(glob.glob("data/raw/rmdb/*.rdat")):
    rid0 = os.path.basename(f)[:-5]
    try:
        r = RDATFile(); r.load(open(f))
    except Exception as e:
        log.append((rid0, f"parse_error:{type(e).__name__}")); continue
    g = r.annotations or {}
    probe = (g.get("modifier") or ["?"])[0]; temp = (g.get("temperature") or [""])[0]
    for ci, c in enumerate(r.constructs.values()):
        seq = c.sequence; ref = c.structure if valid_ref(c.structure, len(seq)) else ""
        wt = [d for d in c.data if (d.annotations.get("mutation") or ["WT"])[0].upper() in ("WT", "")
              and not d.annotations.get("ligand") ]
        for di, d in enumerate(wt[:2]):   # up to two WT replicates
            p = (d.annotations.get("modifier") or [probe])[0]
            if p.lower() in ("nomod", "none", "?"): continue
            dtype = (d.annotations.get("datatype") or ["REACTIVITY"])[0].upper()
            if not dtype.startswith("REACTIVITY"):      # READS etc. are counts, not reactivities
                log.append((rid0, f"datatype_{dtype.lower()}")); continue
            vals = np.array(d.values, float)
            sp = np.array(d.seqpos if d.seqpos else c.seqpos, int)
            if len(sp) != len(vals): log.append((rid0, "seqpos_len_mismatch")); continue
            # signal-to-noise: source flag if present, else mean(value)/mean(error)
            snr_flag = (d.annotations.get("signal_to_noise") or [""])[0].lower()
            if snr_flag.startswith("weak"): log.append((rid0, "snr_weak")); continue
            if d.errors and len(d.errors) == len(vals):
                e = np.array(d.errors, float); snr = np.nanmean(np.clip(vals, 0, None)) / max(np.nanmean(np.abs(e)), 1e-9)
                if snr < MIN_SNR: log.append((rid0, f"snr<{MIN_SNR}")); continue
            add(f"{rid0}_{ci}", seq, p, temp, "RMDB", sp - c.offset, vals, ref, rep=di + 1)

# ---- MERGE-RNA natural DMS-MaP (highest DMS conc minus 0 mM background)
MR = "tools/MERGE-RNA/data/structured_rnas"
for sysdir in sorted(glob.glob(f"{MR}/*")):
    name = os.path.basename(sysdir); ref = ""
    rf = glob.glob(f"{sysdir}/reference/*")
    if rf:
        L = open(rf[0]).read().split()
        ref = next((x for x in L if set(x) <= set("().")), "")
    for rep in (1, 2):
        hi = sorted(glob.glob(f"{sysdir}/Rep{rep}_*mM"), key=lambda p: int(re.findall(r"_(\d+)mM", p)[0]))
        if len(hi) < 2: continue
        d1 = load_rc(f"{hi[-1]}/LAST_MD_sorted.rc", "hi"); d0 = load_rc(f"{hi[0]}/LAST_MD_sorted.rc", "bg")
        seq = "".join(d1.ref_nt).replace("T", "U")
        cov_ok = (d1.total_count >= CFG["min_coverage"]) & (d0.total_count >= CFG["min_coverage"])
        v = np.where(cov_ok, d1.mut_rate - d0.mut_rate, np.nan)
        add(f"MERGE_{name}", seq, "DMS", "37C", "MERGE-RNA", d1.pos.values, v,
            ref if valid_ref(ref, len(seq)) else "", rep=rep)

os.makedirs("data/processed", exist_ok=True)
df = pd.DataFrame(rows, columns=["rna_id", "profile", "pos", "nt", "probe", "reactivity"])
df.to_csv("data/processed/reactivity.tsv.gz", sep="\t", index=False)
R = pd.DataFrame(rnas); R.to_csv("data/processed/rnas.tsv", sep="\t", index=False)
pd.DataFrame(log, columns=["id", "reason"]).to_csv("results/qc_dropped.tsv", sep="\t", index=False)
print(f"profiles={df.profile.nunique()} rnas={len(R)} positions={len(df)} dropped={len(log)}")
print(R.groupby("source").length.describe())
