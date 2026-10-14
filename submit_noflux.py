"""Submit the Figure 10 runs with the no-perturbation-flux buoyancy condition.

Same runs, parameters and resources as the March set (run_all.py), but using
full_code_diagnostics_noflux.jl and writing to OUTPUT_DIR. The SI-unstable runs
(Ri_inv = 1.05) are not included because they are not in the figure.

Usage: python submit_noflux.py [--dry-run] [SIM ...]
  --dry-run   write and print the job scripts but do not qsub them
  SIM ...     submit only these parameter sets (default: all of ALL_SIMS)

Progress is written live to logs/noflux/SIM.log (tail -f it); the PBS .out/.err files
only appear when the job ends.

Precompile the Julia packages once before submitting, so the jobs do not all compile
into the shared depot at the same time:
  julia --pkgimages=no --project=. -e 'using Pkg; Pkg.precompile()'
"""
import subprocess
import sys
from pathlib import Path

ALL_SIMS = [
    "S05R065delta08", "S05R075delta08", "S05R085delta08", "S05R095delta08",   # Sim1-4
    "S10R065delta08", "S10R075delta08", "S10R085delta08", "S10R095delta08",   # Sim6-9
    "S15R075delta08", "S15R085delta08", "S15R095delta08",                     # Sim11-13
    "S10R080delta08", "S10R080delta06", "S10R080delta04", "S10R080delta02",   # Sim15-18
]

OUTPUT_DIR = "/glade/derecho/scratch/wenegrat/BottomBoundaryLayer/output/"
JULIA_FILE = "full_code_diagnostics_noflux.jl"
JOB_DIR = Path("jobs_noflux")
LOG_DIR = Path("logs/noflux")

PBS_SCRIPT = """#!/bin/bash -l
#PBS -N {name}_nf
#PBS -A UMCP0023
#PBS -l select=1:ncpus=32:mpiprocs=2:mem=160gb:ngpus=1:gpu_type=v100
#PBS -l walltime=8:00:00
#PBS -q casper
#PBS -o {log_dir}/{name}.out
#PBS -e {log_dir}/{name}.err
#PBS -m abe

module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.10.2

cd {repo}
# Log straight to glade (PBS only copies the -o/-e files at job end). `script` runs julia
# on a pseudo-terminal so its output is line-buffered; -f flushes each write, -e keeps
# julia's exit status.
script -qfec "julia --pkgimages=no --project=. {julia_file} --path {output_dir} --suffix {name}" {log_dir}/{name}.log
"""


def main(argv):
    dry_run = "--dry-run" in argv
    sims = [a for a in argv if a != "--dry-run"] or ALL_SIMS
    unknown = set(sims) - set(ALL_SIMS)
    if unknown:
        sys.exit(f"not in ALL_SIMS: {sorted(unknown)}")

    repo = Path(__file__).resolve().parent
    for d in (repo / JOB_DIR, repo / LOG_DIR, Path(OUTPUT_DIR)):
        d.mkdir(parents=True, exist_ok=True)

    for name in sims:
        job = repo / JOB_DIR / f"{name}.sh"
        job.write_text(PBS_SCRIPT.format(name=name, log_dir=repo / LOG_DIR, repo=repo,
                                         julia_file=JULIA_FILE, output_dir=OUTPUT_DIR))
        if dry_run:
            print(f"--- {job}\n{job.read_text()}")
        else:
            out = subprocess.run(["qsub", str(job)], capture_output=True, text=True, check=True)
            print(f"{name}: {out.stdout.strip()}")


if __name__ == "__main__":
    main(sys.argv[1:])
