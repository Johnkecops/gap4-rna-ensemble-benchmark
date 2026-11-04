"""Select RMDB entries for natural/structured RNAs: drop Eterna category and
high-throughput design libraries. Prints '<rmdb_id> <rdat_url>'."""
import re, sys, pathlib
SKIP = re.compile(r"^(ETERNA|15K|OPENKNOT|PK50|PK90|RIBONANZA|RYOS|FSEATTR|SYN|LIB)", re.I)
# >20 MB high-throughput libraries, removed on user decision at 55% disk (2026-11-04)
EXCL_F = pathlib.Path(sys.argv[1]).parent.parent / "excluded_large_libraries.txt"
EXCL = set(EXCL_F.read_text().split()) if EXCL_F.exists() else set()
for f in sorted(pathlib.Path(sys.argv[1]).glob("*.md")):
    t = f.read_text(errors="ignore")
    rid = re.search(r'rmdb_id:\s*"([^"]+)"', t); url = re.search(r"^rdat:\s*(\S+)", t, re.M)
    cat = re.search(r'category:\s*"([^"]+)"', t)
    if not (rid and url) or (cat and cat.group(1) == "Eterna") or "data-eterna" in url.group(1):
        continue
    if SKIP.match(rid.group(1)) or rid.group(1) in EXCL:
        continue
    print(rid.group(1), url.group(1))
