"""Energy-model consistency check: infinite-beam linear tools must reproduce ViennaRNA."""
import sys, os, random, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import foldtools as ft

random.seed(1)
SEQS = ["".join(random.choice("ACGU") for _ in range(100)) for _ in range(5)]


def test_linearpartition_inf_beam_matches_vienna():
    for s in SEQS:
        _, pv, _ = ft.vienna(s, k=0)
        _, pl = ft.linearpartition(s, 0)
        assert np.max(np.abs(pv - pl)) < 0.01, np.max(np.abs(pv - pl))


def test_linearsampling_exact_matches_vienna():
    k = 2000
    for s in SEQS:
        _, pv, _ = ft.vienna(s, k=0)
        ps = ft.punp_from_samples(ft.linearsampling(s, 0, k), len(s))
        # 5 binomial SE at p=0.5 with k samples ~ 0.056
        assert np.max(np.abs(pv - ps)) < 5 * np.sqrt(0.25 / k), np.max(np.abs(pv - ps))
