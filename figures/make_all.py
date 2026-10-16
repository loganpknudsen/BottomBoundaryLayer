"""Generate every manuscript figure (figures/NN_*.py) into figures/output/.

Each script runs in its own process, so one failure (e.g. a simulation that has not finished)
does not stop the rest; a summary at the end lists what failed, with the last lines of its error.
HDF5 file locking is switched off so the scripts can safely read files that running simulations
are still writing. Data locations come from paths.py (BBL_SIM_DIR, BBL_LT_DIR, BBL_FIG_OUT).

Usage:
    python figures/make_all.py              # all figures
    python figures/make_all.py 10 14 17     # only these
    python figures/make_all.py --residual   # also the budget-residual versions of Figs 7, 14, 16
"""
import os
import subprocess
import sys
import time
from pathlib import Path

FIG_DIR = Path(__file__).resolve().parent
RESIDUAL_FIGS = {7, 14, 16}

args = sys.argv[1:]
with_residual = "--residual" in args
wanted = {int(a) for a in args if a != "--residual"}

scripts = sorted(p for p in FIG_DIR.glob("[0-9][0-9]_*.py"))
if wanted:
    scripts = [p for p in scripts if int(p.name[:2]) in wanted]
jobs = [(p, []) for p in scripts]
if with_residual:
    jobs += [(p, ["--residual"]) for p in scripts if int(p.name[:2]) in RESIDUAL_FIGS]
    jobs.sort(key=lambda j: j[0].name)

env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE", MPLBACKEND="Agg")
failed = []
for script, extra in jobs:
    name = " ".join([script.name, *extra])
    print(f"=== {name}", flush=True)
    t0 = time.time()
    result = subprocess.run([sys.executable, str(script), *extra], env=env, cwd=FIG_DIR,
                            capture_output=True, text=True)
    print(result.stdout, end="")
    if result.returncode == 0:
        print(f"    ok ({time.time() - t0:.0f} s)", flush=True)
    else:
        print(f"    FAILED ({time.time() - t0:.0f} s)", flush=True)
        failed.append((name, result.stderr.strip().splitlines()[-3:]))

print(f"\n{len(jobs) - len(failed)} of {len(jobs)} succeeded")
for name, err in failed:
    print(f"  {name}:")
    for line in err:
        print(f"      {line}")
sys.exit(1 if failed else 0)
