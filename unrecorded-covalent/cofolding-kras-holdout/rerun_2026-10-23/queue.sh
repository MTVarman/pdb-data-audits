#!/bin/bash
# Re-run of the arms whose per-sample data was lost on 22 Oct 2026.
# Order is by how much of the paper depends on each arm.
cd /content
export LAYERNORM_TYPE=torch
export DEST=${DEST:-/content/drive/MyDrive/kras_cofold_rerun_2026-10-23}
PX=/content/px/bin/protenix
BZ=/content/bz/bin/boltz
log(){ echo "[$(date -u +%H:%M:%S)Z] $*" | tee -a /content/queue.log; }
PXFLAGS="-s 42 -e 5 -p 200 -c 10 --use_default_params true"
BZFLAGS="--diffusion_samples 5 --seed 42 --no_kernels --output_format mmcif --use_msa_server"

run_px () { # run_px MODEL INPUT_JSON OUT_RAW_DIR
  local o=$3
  [ "$(find $o -name '*.cif' 2>/dev/null | wc -l)" -gt 0 ] && { log "skip $(basename $o)"; return; }
  mkdir -p $o; log "px:$1 $(basename $o)"
  $PX pred -i $2 -o $o -n $1 $PXFLAGS >> /content/s_px.log 2>&1
  log "   -> $(find $o -name '*.cif' | wc -l) cifs"
}
run_bz () { # run_bz YAML OUT_RAW_DIR
  local o=$2
  [ "$(find $o -name '*.cif' 2>/dev/null | wc -l)" -gt 0 ] && { log "skip $(basename $o)"; return; }
  mkdir -p $o; log "bz $(basename $o)"
  $BZ predict $1 --out_dir $o $BZFLAGS >> /content/s_bz.log 2>&1
  local n=$(find $o -name '*.cif' | wc -l)
  if [ "$n" -eq 0 ]; then log "   -> FAILED: $(grep -E 'Error|error|not found' /content/s_bz.log | tail -1 | cut -c1-160)"; else log "   -> $n cifs"; fi
}

log "=== PREPARE ==="
python kras_holdout.py prepare work_px --manifest run_manifest.csv 2>&1 | tee -a /content/queue.log | tail -3
python kras_holdout.py prepare work_control --manifest control_manifest.csv --tier D 2>&1 | tee -a /content/queue.log | tail -3
python protenix_arm.py prepare work_px 2>&1 | tail -2 | tee -a /content/queue.log
python mk_control_px.py work_control 2>&1 | tail -2 | tee -a /content/queue.log
python check_px_bonds.py work_px 2>&1 | tee -a /content/queue.log | tail -3 || { log "BOND CHECK FAILED - stopping"; exit 1; }
test $(ls work_px/inputs/protenix/*.json | wc -l) -eq 42 || { log "expected 42 protenix inputs"; exit 1; }
test $(ls work_control/inputs/protenix/*.json | wc -l) -eq 40 || { log "expected 40 control inputs"; exit 1; }
mkdir -p work_px_20250630/structures && cp work_px/manifest.json work_px_20250630/ && cp work_px/structures/* work_px_20250630/structures/
cp -r work_px/inputs work_px_20250630/inputs
mkdir -p "$DEST" && cp -f work_px/manifest.json "$DEST/work_px_manifest.json" && cp -f work_control/manifest.json "$DEST/work_control_manifest.json"
log "prepare done"

log "=== STAGE 1a: protenix-20250630, Tier B + C, both arms (30 jobs) ==="
for e in 9HMR 9I7Y 9N44 9NZM 9NZN 9PZF 9PZY 8S8C 8V39 8V3A 9DMM 9E9H 9E9I 9KPM 9KPN; do
  for a in smiles bond; do run_px protenix_base_20250630_v1.0.0 work_px/inputs/protenix/${e}_${a}.json work_px_20250630/out_raw/protenix/${e}_${a}; done
done
python collect.py work_px_20250630 protenix | tail -1 | tee -a /content/queue.log
bash persist.sh S1a work_px_20250630 | tee -a /content/queue.log
log "STAGE 1a COMPLETE"

log "=== STAGE 1b: protenix-20250630, Tier A, both arms (12 jobs) ==="
for e in 9BAI 9BAJ 9BAK 8T4V 9GBJ 9QQ0; do
  for a in smiles bond; do run_px protenix_base_20250630_v1.0.0 work_px/inputs/protenix/${e}_${a}.json work_px_20250630/out_raw/protenix/${e}_${a}; done
done
python collect.py work_px_20250630 protenix | tail -1 | tee -a /content/queue.log
bash persist.sh S1b work_px_20250630 | tee -a /content/queue.log
log "STAGE 1b COMPLETE"

log "=== STAGE 2: boltz-2, all 21 entries, both arms (42 jobs; CCD-missing bond arms fail by design) ==="
for d in work_px/inputs/boltz/*/; do
  n=$(basename $d); y=$(ls $d*.yaml 2>/dev/null | head -1); [ -z "$y" ] && continue
  run_bz $y work_px/out_raw/boltz/$n
done
python collect.py work_px boltz | tail -1 | tee -a /content/queue.log
bash persist.sh S2 work_px | tee -a /content/queue.log
log "STAGE 2 COMPLETE"

log "=== STAGE 3: boltz-2, 40 controls, smiles only ==="
for d in work_control/inputs/boltz/*_smiles/; do
  n=$(basename $d); y=$(ls $d*.yaml 2>/dev/null | head -1); [ -z "$y" ] && continue
  run_bz $y work_control/out_raw/boltz/$n
done
python collect.py work_control boltz | tail -1 | tee -a /content/queue.log
bash persist.sh S3 work_control | tee -a /content/queue.log
log "STAGE 3 COMPLETE"

log "=== STAGE 4: protenix v1.0.0, 40 controls, smiles only ==="
for f in work_control/inputs/protenix/*_smiles.json; do
  n=$(basename $f .json); run_px protenix_base_default_v1.0.0 $f work_control/out_raw/protenix/$n
done
python collect.py work_control protenix | tail -1 | tee -a /content/queue.log
bash persist.sh S4 work_control | tee -a /content/queue.log
log "STAGE 4 COMPLETE"
log "ALL STAGES COMPLETE"
