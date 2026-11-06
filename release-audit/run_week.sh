#!/usr/bin/env bash
# Audit the most recent PDB release (Wednesdays) and file it under weeks/.
# Usage: ./run_week.sh [YYYY-MM-DD]
set -euo pipefail
cd "$(dirname "$0")"
DATE="${1:-$(python3 -c 'import datetime as d; t=d.date.today(); print((t-d.timedelta(days=(t.weekday()-2)%7)).isoformat())')}"
WORK="${TMPDIR:-/tmp}/pdb_release_$DATE"
python3 -I release_audit.py "$DATE" "$WORK"
mkdir -p weeks
cp "$WORK/${DATE}_findings.csv" "$WORK/${DATE}_brief.md" weeks/
echo "filed weeks/${DATE}_findings.csv"
