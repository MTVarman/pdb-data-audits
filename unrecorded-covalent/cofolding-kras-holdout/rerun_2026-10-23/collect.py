"""Flatten raw model output into out/<model>/<entry>_<arm>/sample_i.cif.

Copies from out_raw/ (outside out/) so a recursive glob over out/ sees each structure once —
the 21 Oct Boltz Tier C double count came from copying model_i.cif next to the
nested boltz_results/ tree and then rglob-ing both.
"""
import shutil, sys, pathlib
work = pathlib.Path(sys.argv[1]); model = sys.argv[2]
raw = work/"out_raw"/model
if not raw.exists():
    print("no out_raw at", raw); sys.exit(0)
n = 0
for job in sorted(p for p in raw.iterdir() if p.is_dir()):
    dest = work/"out"/model/job.name
    cifs = sorted(job.rglob("*.cif"))
    if not cifs:
        continue
    dest.mkdir(parents=True, exist_ok=True)
    for i, c in enumerate(cifs):
        shutil.copy(c, dest/f"sample_{i}.cif"); n += 1
    print(f"  {job.name}: {len(cifs)}")
print(f"collected {n} structures into {work}/out/{model}")
