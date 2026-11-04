"""Checklist for RNA-structure-Manuscript-Arli.docx (stand-in for the /goal loop). Exit 1 if any check fails."""
import re, sys, os, csv
from docx import Document
from docx.oxml.ns import qn

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "..", "RNA-structure-Manuscript-Arli.docx")
d = Document(F)
M = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"
res = []
def check(name, ok, info=""): res.append((ok, name, info))

# walk top-level body paragraphs with their section
paras, sec = [], "front"
for p in d.paragraphs:
    t = p.text.strip(); st = p.style.name
    if st.startswith("Heading"):
        sec = t; paras.append((sec, "H", t, p)); continue
    has_math = p._p.find(f".//{M}oMathPara") is not None
    kind = "EQ" if has_math and not re.sub(r"\s", "", t) else ("CAP" if st in ("Image Caption", "Table Caption", "Captioned Figure") or t.startswith(("Fig. ", "Table ")) else "P")
    if t or has_math: paras.append((sec, kind, t, p))

heads = [t for s, k, t, p in paras if k == "H"]
need = ["Highlights", "Abstract", "1. Introduction", "2. Materials and methods", "3. Results", "4. Discussion", "5. Conclusions",
        "CRediT authorship contribution statement", "Declaration of competing interest",
        "Declaration of generative AI and AI-assisted technologies in the writing process", "Data and code availability",
        "Funding", "Acknowledgements", "References"]
idx = [heads.index(h) if h in heads else -1 for h in need]
check("required sections present and ordered", -1 not in idx and idx == sorted(idx), str([h for h, i in zip(need, idx) if i < 0]))

hl = [t for s, k, t, p in paras if s == "Highlights" and k == "P"]
check("highlights: 3 to 5 bullets, each <= 85 characters", 3 <= len(hl) <= 5 and all(len(x) <= 85 for x in hl), str([len(x) for x in hl]))
ab = " ".join(t for s, k, t, p in paras if s == "Abstract" and k == "P" and not t.startswith("Keywords"))
check("abstract <= 250 words", len(ab.split()) <= 250, f"{len(ab.split())} words")
kw = next((t for s, k, t, p in paras if t.startswith("Keywords")), "")
nkw = len(kw.split(":", 1)[-1].split(";")) if kw else 0
check("keywords 1 to 7", 1 <= nkw <= 7, f"{nkw}")

CITE = re.compile(r"\[(\d+(?:,\d+)*)\]")
check("no citation in abstract or keywords", not CITE.search(ab + kw))
concl = " ".join(t for s, k, t, p in paras if s == "5. Conclusions")
check("no citation in conclusions", not CITE.search(concl))

body_secs = lambda s: re.match(r"^[1-4]\.", s) is not None
missing = [t[:70] for s, k, t, p in paras if body_secs(s) and k in ("P", "CAP") and not CITE.search(t)]
check("every body paragraph and caption (sections 1 to 4) cites at least one reference", not missing, "; ".join(missing))

# equations and pseudo-code: explanatory paragraph right before each equation must cite
eq_ok, bad = True, []
for i, (s, k, t, p) in enumerate(paras):
    if k == "EQ":
        prev = next((paras[j][2] for j in range(i - 1, -1, -1) if paras[j][1] == "P"), "")
        nxt = next((paras[j][2] for j in range(i + 1, len(paras)) if paras[j][1] in ("P", "H")), "")
        if not (CITE.search(prev) or CITE.search(nxt)): eq_ok = False; bad.append(prev[:50])
neq = sum(1 for x in paras if x[1] == "EQ")
check("display equations are Office Math (>= 12) and each is explained with a citation", neq >= 12 and eq_ok, f"{neq} equations; {bad}")
ninline = len(d.element.body.findall(f".//{M}oMath")) - neq
check("inline symbols set as Office Math", ninline >= 20, f"{ninline} inline math objects")
uni = all(re.search(r"[ρτΣ∑≈∈ΔΩβμ√]", "".join(x.text or "" for x in d.element.body.iter(f"{M}t"))) for _ in [0])
check("equations contain Unicode math characters", uni)

boxes = [t for t in d.tables if t.cell(0, 0).paragraphs[0].text.startswith("Algorithm")]
bordered = all(t._tbl.tblPr.find(qn("w:tblBorders")) is not None for t in boxes)
check("pseudo-code in 4 bordered boxes", len(boxes) == 4 and bordered, f"{len(boxes)} boxes")
alg_cited = all(any(re.search(rf"Algorithm {n}", t) and CITE.search(t) for s, k, t, p in paras if k == "P") for n in (1, 2, 3, 4))
check("every algorithm is introduced in a cited paragraph", alg_cited)

# references
refp = [t for s, k, t, p in paras if s == "References" and k == "P"]
nref = len(refp)
cited = set()
for s, k, t, p in paras:
    if s not in ("References",):
        for g in CITE.findall(t): cited |= {int(x) for x in g.split(",")}
for tb in d.tables:
    for row in tb.rows:
        for c in row.cells:
            for g in CITE.findall(c.text): cited |= {int(x) for x in g.split(",")}
check("every reference cited in the body, every citation resolves", cited == set(range(1, nref + 1)), f"refs {nref}, uncited {sorted(set(range(1, nref + 1)) - cited)}, dangling {sorted(cited - set(range(1, nref + 1)))}")
val = {r["pmid"]: r for r in csv.DictReader(open(os.path.join(HERE, "reference_validation.tsv")), delimiter="\t")}
okref = []
for t in refp:
    pm = re.search(r"PMID: (\d+)", t); doi = re.search(r"https://doi.org/(\S+?)\. PMID", t)
    okref.append(bool(pm and doi and pm.group(1) in val and val[pm.group(1)]["doi"] == doi.group(1) and val[pm.group(1)]["pubmed_doi_match"] == "yes"))
check("every reference has a PubMed-validated PMID and matching DOI", all(okref) and nref == len(val), f"{sum(okref)}/{nref} validated")

# numbers cross-checked against results/*.tsv
import pandas as pd
R = os.path.join(HERE, "..", "results")
full = "\n".join(p.text for p in d.paragraphs) + "\n".join(c.text for t in d.tables for r in t.rows for c in r.cells)
L = pd.read_csv(f"{R}/leaderboard.tsv", sep="\t").set_index(["stratum", "tool"])
B = pd.read_csv(f"{R}/bias_summary.tsv", sep="\t").set_index(["tool", "stratum"])
nums = {
 "rho <200 vienna": f"{L.loc[('<200', 'vienna_b-1'), 'spearman']:.3f}",
 "rho 200-1000 vienna": f"{L.loc[('200-1000', 'vienna_b-1'), 'spearman']:.3f}",
 "LP b100 MAD long": f"{B.loc[('linearpartition_b100', '1000-3000'), 'mad_unp_mean']:.3f}",
 "LS b100 MAD long": f"{B.loc[('linearsampling_b100', '1000-3000'), 'mad_unp_mean']:.3f}",
 "positive control": open(f"{R}/positive_control_eternabench.txt").read().split("pearson_ours=")[1][:6],
 "RNAs": "439", "profiles": "614", "positions": "75,568", "max length": "2,216",
}
miss = [k for k, v in nums.items() if v not in full]
check("key numbers match results tables", not miss, str(miss))

# hygiene
hidden = re.findall(r"[​‌‍⁠﻿­]", full)
check("no zero-width or hidden watermark characters", not hidden, f"{len(hidden)} found")
prose = "\n".join(t for s, k, t, p in paras if s not in ("References", "CRediT authorship contribution statement"))
check("no em dashes in prose", "—" not in prose)
cp = d.core_properties
check("document properties clean", cp.author == "Arli Aditya Parikesit" and not cp.comments and "pandoc" not in (cp.last_modified_by or "").lower(), f"{cp.author!r} {cp.last_modified_by!r}")
check("file saved under required name", os.path.basename(F) == "RNA-structure-Manuscript-Arli.docx" and os.path.exists(F))

for ok, n, info in res: print(("PASS " if ok else "FAIL ") + n + (f"  [{info}]" if info and not ok else (f"  ({info})" if info else "")))
print(f"\n{sum(o for o, *_ in res)}/{len(res)} checks passed")
sys.exit(0 if all(o for o, *_ in res) else 1)
