#!/bin/bash
# Stage A: acquire data. Writes data/raw/MANIFEST.tsv (file, sha256, date, source).
# Cloud-session constraints (checked 2026-11-04): GitHub release listing/API, Kaggle,
# Zenodo, RASP and bpRNA hosts are blocked by egress policy. RMDB RDATs are fetched
# one by one from the release URL written in each _entries/*.md of the site repo.
set -euo pipefail
cd "$(dirname "$0")/.."
RAW=data/raw; mkdir -p $RAW/rmdb
[ -d $RAW/EternaBench ] || git clone -q --depth 1 https://github.com/eternagame/EternaBench $RAW/EternaBench
if [ ! -d $RAW/rmdb_repo/_entries ]; then
  git clone -q --depth 1 --filter=blob:none --no-checkout https://github.com/DasLab/rmdb.github.io $RAW/rmdb_repo
  git -C $RAW/rmdb_repo checkout -q HEAD -- _entries
fi
# entry selection (natural/structured RNAs, no Eterna/library sets) done in python
python scripts/select_rmdb.py $RAW/rmdb_repo/_entries > $RAW/rmdb/selected_urls.txt
while read -r id url; do
  f=$RAW/rmdb/$id.rdat
  [ -s "$f" ] || curl -sfL -m 120 -o "$f" "$url" || { echo "FAIL $id" >&2; rm -f "$f"; }
  scripts/diskcheck.sh >/dev/null || { [ $? -eq 2 ] && { echo "disk >=95%, abort"; exit 2; }; }
done < $RAW/rmdb/selected_urls.txt
# MERGE-RNA DMS data (natural RNAs) live inside the cloned tool repo
ln -sfn ../../tools/MERGE-RNA/data $RAW/merge_rna
d=$(date -u +%F)
{ echo -e "file\tsha256\tdownload_date\tsource"
  find $RAW/rmdb -name '*.rdat' | sort | while read f; do echo -e "$f\t$(sha256sum $f|cut -c1-64)\t$d\tRMDB release asset"; done
  for f in $RAW/EternaBench/data/*.zip $RAW/EternaBench/data/ExternalData/*.zip; do echo -e "$f\t$(sha256sum $f|cut -c1-64)\t$d\tEternaBench@$(git -C $RAW/EternaBench rev-parse --short HEAD)"; done
  find -L $RAW/merge_rna/structured_rnas $RAW/merge_rna/cspA -name '*.rc' | sort | while read f; do echo -e "$f\t$(sha256sum $f|cut -c1-64)\t$d\tMERGE-RNA@$(git -C tools/MERGE-RNA rev-parse --short HEAD)"; done
} > $RAW/MANIFEST.tsv
echo "manifest rows: $(($(wc -l < $RAW/MANIFEST.tsv)-1))"
