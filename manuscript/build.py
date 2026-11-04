"""Build RNA-structure-Manuscript-Arli.docx from manuscript.md.
1) number [@key] citations by first appearance (Elsevier numbered style) and write the reference list
2) pandoc -> docx (LaTeX math becomes native Office Math / OMML)
3) python-docx: Times New Roman 12 pt, double spacing, A4, 2.5 cm margins, line and page numbers,
   pseudo-code blocks in bordered boxes, compact tables."""
import re, subprocess, sys, os
from copy import deepcopy
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_LINE_SPACING, WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from refs import REFS, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "RNA-structure-Manuscript-Arli.docx")
md = open(os.path.join(HERE, "manuscript.md"), encoding="utf-8").read()

# ---- 1. citations
order = []
for grp in re.findall(r"\[(@[^\]]+)\]", md):
    for k in re.findall(r"@([\w]+)", grp):
        if k not in REFS: sys.exit(f"unknown citation key {k}")
        if k not in order: order.append(k)
num = {k: i + 1 for i, k in enumerate(order)}
md = re.sub(r"\[(@[^\]]+)\]", lambda m: "[" + ",".join(str(n) for n in sorted({num[k] for k in re.findall(r"@(\w+)", m.group(1))})) + "]", md)
unused = set(REFS) - set(order)
if unused: sys.exit(f"references never cited: {unused}")
reflist = "\n\n".join(f"[{num[k]}] {fmt(k)}" for k in order)
md = md.replace("::: refs\n:::", "::: refs\n" + reflist.replace("[", "\\[").replace("]", "\\]") + "\n:::")
open(os.path.join(HERE, "manuscript_numbered.md"), "w", encoding="utf-8").write(md)

# ---- 2. pandoc
tmp = os.path.join(HERE, "_pandoc.docx")
subprocess.run(["pandoc", "manuscript_numbered.md", "-f", "markdown+tex_math_dollars+fenced_divs", "-o", tmp,
                "--resource-path", HERE], cwd=HERE, check=True)

# ---- 3. formatting
d = Document(tmp)
for name in ("Normal", "Body Text", "First Paragraph", "Compact", "Image Caption", "Table Caption", "Caption"):
    try:
        st = d.styles[name]
    except KeyError:
        continue
    st.font.name = "Times New Roman"; st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
for name in ("Title", "Heading 1", "Heading 2", "Heading 3"):
    try:
        st = d.styles[name]
    except KeyError:
        continue
    st.font.name = "Times New Roman"; st.font.color.rgb = RGBColor(0, 0, 0); st.font.bold = True
    st.font.size = Pt({"Title": 14, "Heading 1": 14, "Heading 2": 12, "Heading 3": 12}[name])
    st.font.italic = name == "Heading 3"
    st.element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman"); st.element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")

body = d.element.body
for par in d.paragraphs:   # direct run formatting beats theme fonts and colours in pandoc styles
    if par.style.name in ("Title", "Heading 1", "Heading 2", "Heading 3"):
        for r in par.runs:
            r.font.name = "Times New Roman"; r.font.color.rgb = RGBColor(0, 0, 0); r.font.bold = True
            r.font.size = Pt(14 if par.style.name in ("Title", "Heading 1") else 12)
            rf = r._r.get_or_add_rPr().get_or_add_rFonts()
            for a in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"): rf.set(qn(a), "Times New Roman")
            for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
                if rf.get(qn(a)) is not None: del rf.attrib[qn(a)]


def make_box(par):
    """replace a pseudo-code paragraph by a 1x1 bordered table holding monospaced lines"""
    text = "".join(t.text or "" for t in par._p.iter(qn("w:t")))
    lines = []
    for r in par._p.iter(qn("w:r")):
        for ch in r:
            if ch.tag == qn("w:t"): lines.append(("t", ch.text or ""))
            elif ch.tag == qn("w:br"): lines.append(("br", ""))
    tbl = d.add_table(rows=1, cols=1)
    cell = tbl.cell(0, 0)
    tblPr = tbl._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for side in ("top", "left", "bottom", "right"):
        e = OxmlElement(f"w:{side}"); e.set(qn("w:val"), "single"); e.set(qn("w:sz"), "8"); e.set(qn("w:color"), "000000")
        borders.append(e)
    tblPr.append(borders)
    shd = OxmlElement("w:shd"); shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), "F2F2F2")
    cell._tc.get_or_add_tcPr().append(shd)
    content = "".join(x if k == "t" else "\n" for k, x in lines)
    first = True
    for ln in content.split("\n"):
        p = cell.paragraphs[0] if first else cell.add_paragraph()
        r = p.add_run(ln); r.font.name = "Courier New"; r.font.size = Pt(8.5)
        r._r.get_or_add_rPr().get_or_add_rFonts().set(qn("w:eastAsia"), "Courier New")
        if first: r.bold = True
        p.paragraph_format.space_after = Pt(0); p.paragraph_format.line_spacing = 1.0
        first = False
    par._p.addprevious(tbl._tbl)
    par._p.getparent().remove(par._p)


for par in list(d.paragraphs):
    if par.style.name == "Source Code" or par._p.find(qn("w:pPr")) is not None and "SourceCode" in (par._p.pPr.xml if par._p.pPr is not None else ""):
        make_box(par)

for par in d.paragraphs:
    pf = par.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    pf.space_after = Pt(0) if par.style.name.startswith("Heading") is False else Pt(6)
    if par.style.name in ("Body Text", "First Paragraph"):
        pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
in_refs = False
for par in d.paragraphs:   # reference list: left aligned
    if par.style.name.startswith("Heading"): in_refs = par.text.strip() == "References"
    elif in_refs: par.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
for tbl in d.tables:
    for row in tbl.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                p.paragraph_format.line_spacing = 1.0; p.paragraph_format.space_after = Pt(0)
                for r in p.runs:
                    if r.font.name != "Courier New": r.font.size = Pt(10); r.font.name = "Times New Roman"

# page setup, line numbers, page numbers
for s in d.sections:
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Cm(2.5)
    ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "283")
    s._sectPr.append(ln)
    fp = s.footer.paragraphs[0]; fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = fp.add_run()
    for tag, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), tag)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt
        run._r.append(e)

cp = d.core_properties
cp.author = "Arli Aditya Parikesit"; cp.last_modified_by = "Arli Aditya Parikesit"
cp.title = "Benchmarking RNA secondary structure ensemble generators against chemical mapping of natural RNAs"
cp.comments = ""; cp.keywords = "RNA secondary structure; partition function; chemical probing"; cp.subject = ""; cp.category = ""
d.save(OUT); os.remove(tmp)
print("saved", OUT, "| references:", len(order), "| boxes:", sum(1 for t in d.tables if t.cell(0, 0).paragraphs[0].runs and t.cell(0, 0).paragraphs[0].runs[0].font.name == "Courier New"))

# replace pandoc's template app.xml (stale word counts, fixed application string) with a minimal part
import zipfile, shutil
APP = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Properties xmlns="http://schemas.openxmlformats.org/'
       'officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
       '<DocSecurity>0</DocSecurity><SharedDoc>false</SharedDoc></Properties>')
tmpz = OUT + ".tmp"
with zipfile.ZipFile(OUT) as zin, zipfile.ZipFile(tmpz, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        zout.writestr(item, APP if item.filename == "docProps/app.xml" else zin.read(item.filename))
shutil.move(tmpz, OUT)
