#!/bin/bash
# Submit the linear-theory parameter grids (corrected Eq 2.27 system) as Casper CPU jobs.
# Usage (from the repo root): bash linear_theory/submit_grids.sh [sgamma] [freq]
# Output: $OUT/<grid>.nc; live progress in logs/linear_theory/<grid>.log

OUT=/glade/derecho/scratch/$USER/BottomBoundaryLayer/output/linear_theory
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT" "$REPO/logs/linear_theory"
GRIDS=${@:-sgamma freq}

for g in $GRIDS; do
qsub <<EOF
#!/bin/bash -l
#PBS -N lt_$g
#PBS -A UMCP0023
#PBS -q casper
#PBS -l select=1:ncpus=36:mem=100GB
#PBS -l walltime=03:00:00
#PBS -o $REPO/logs/linear_theory/$g.pbs.out
#PBS -e $REPO/logs/linear_theory/$g.pbs.err

module --force purge
module load ncarenv conda
conda activate npl
cd $REPO/linear_theory
python -u stability_grid.py $g 0 1 $OUT 36 > $REPO/logs/linear_theory/$g.log 2>&1 &&
python -u merge_grid.py $g $OUT >> $REPO/logs/linear_theory/$g.log 2>&1
EOF
done
