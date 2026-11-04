#!/bin/bash
# prints used% of session allowance; exit 2 if >=95, 1 if >=50
read used avail < <(df -B1 --output=used,avail / | tail -1)
pct=$(( 100*used/(used+avail) ))
echo "disk_used_pct=$pct used=$((used/1073741824))G avail=$((avail/1073741824))G"
[ $pct -ge 95 ] && exit 2; [ $pct -ge 50 ] && exit 1; exit 0
