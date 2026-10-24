#!/bin/bash
# persist.sh STAGE WORKDIR  — score, pose-RMSD, then copy everything that matters to Drive.
STAGE=$1; WORK=$2
cd /content
DEST=${DEST:-/content/drive/MyDrive/kras_cofold_rerun_2026-10-23}
mkdir -p "$DEST/$WORK"
python /content/kras_holdout.py score $WORK > $WORK/score_${STAGE}.txt 2>&1 || true
python /content/poserms.py $WORK > $WORK/poserms_${STAGE}.txt 2>&1 || true
cp -f $WORK/results.csv   "$DEST/$WORK/results.csv"   2>/dev/null
cp -f $WORK/pose_rmsd.csv "$DEST/$WORK/pose_rmsd.csv" 2>/dev/null
cp -f $WORK/manifest.json "$DEST/$WORK/manifest.json" 2>/dev/null
cp -f $WORK/score_${STAGE}.txt $WORK/poserms_${STAGE}.txt "$DEST/$WORK/" 2>/dev/null
tar -czf "$DEST/${WORK}_out_${STAGE}.tar.gz" -C /content $WORK/out $WORK/inputs 2>/dev/null
cp -f /content/queue.log /content/s*.log "$DEST/" 2>/dev/null
echo "$(date -u +%FT%TZ) $STAGE $WORK rows=$(( $(wc -l < $WORK/results.csv 2>/dev/null || echo 1) - 1 ))" >> "$DEST/STAGES_COMPLETE.txt"
sync
echo "persisted $STAGE -> $DEST"
