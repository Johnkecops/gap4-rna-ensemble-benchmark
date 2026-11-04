"""Thin wrappers: each returns (bpp dict {(i,j):p} 1-based, punp array, samples list)."""
import subprocess, os, numpy as np, RNA
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LP = f"{ROOT}/tools/LinearPartition/bin/linearpartition_v"
LS = f"{ROOT}/tools/LinearSampling/bin"
GTIME = os.path.join(os.path.dirname(os.path.realpath(__import__("sys").executable)), "time")
LAST_RSS_MB = [0.0]


def timed(cmd, inp):
    """run cmd under GNU time; store child peak RSS (MB) in LAST_RSS_MB[0]"""
    tf = f"/tmp/claude-0/rss_{os.getpid()}"
    r = subprocess.run([GTIME, "-f", "%M", "-o", tf] + cmd, input=inp, text=True, capture_output=True, check=True)
    LAST_RSS_MB[0] = int(open(tf).read().split()[-1]) / 1024; os.remove(tf)
    return r


def punp_from_bpp(bpp, n):
    pu = np.ones(n)
    for (i, j), p in bpp.items():
        pu[i - 1] -= p; pu[j - 1] -= p
    return np.clip(pu, 0, 1)


def vienna(seq, k=1000, cutoff=1e-4):
    md = RNA.md(); md.uniq_ML = 1
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe(); fc.exp_params_rescale(mfe); fc.pf()
    P = fc.bpp(); n = len(seq)
    bpp = {(i, j): P[i][j] for i in range(1, n + 1) for j in range(i + 1, n + 1) if P[i][j] > cutoff}
    samples = list(fc.pbacktrack(k)) if k else []
    return bpp, punp_from_bpp(bpp, n), samples


def linearpartition(seq, beam, cutoff=1e-4, tmp="/tmp/claude-0"):
    out = f"{tmp}/lp_{os.getpid()}.bpp"
    if os.path.exists(out): os.remove(out)
    # args: beam sharpturn verbose bpp_file prefix pf_only cutoff forest mea gamma TK thr TKprefix MEAprefix bpseq shape fasta dangles evaly
    cmd = [LP, str(beam), "0", "0", out, "", "0", str(cutoff), "", "0", "3.0", "0", "0.3", "", "", "0", "", "0", "2", ""]
    timed(cmd, seq + "\n")
    bpp = {}
    for line in open(out):
        p = line.split()
        if len(p) == 3: bpp[(int(p[0]), int(p[1]))] = float(p[2])
    os.remove(out)
    return bpp, punp_from_bpp(bpp, len(seq))


def linearsampling(seq, beam, k=1000):
    exe = f"{LS}/{'exact_' if beam <= 0 else ''}linearsampling_lazysaving"
    # args: beam sharpturn verbose k readforest fasta input_file shape
    r = timed([exe, str(beam), "0", "0", str(k), "", "0", "", ""], seq + "\n")
    return [l.strip() for l in r.stdout.splitlines() if l and l[0] in "(.)" and len(l.strip()) == len(seq)]


def punp_from_samples(samples, n):
    pu = np.zeros(n)
    for s in samples: pu += np.frombuffer(s.encode(), np.uint8) == ord(".")
    return pu / max(len(samples), 1)
