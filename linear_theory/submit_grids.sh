#!/bin/bash
# Submit the linear-theory parameter grids (corrected Eq 2.27 system) as Casper CPU jobs.
# Usage (from the repo root): [REFINE=n] bash linear_theory/submit_grids.sh [sgamma] [freq]
#   REFINE=n divides the grid spacing by n (default 1, the preprint grids); cost ~ n^2.
# Output: $OUT/<grid>.nc, with OUT = output/linear_theory (REFINE=1) or output/linear_theory_r<n>;
# live progress in logs/linear_theory/<grid>_r<n>.log

REFINE=${REFINE:-1}
if [ "$REFINE" = 1 ]; then SUFFIX=""; else SUFFIX="_r$REFINE"; fi
OUT=/glade/derecho/scratch/$USER/BottomBoundaryLayer/output/linear_theory$SUFFIX
REPO=$(cd "$(dirname "$0")/.." && pwd)
LOGS=$REPO/logs/linear_theory
mkdir -p "$OUT" "$LOGS"
GRIDS=${@:-sgamma freq}
HOURS=$(( 3 * REFINE * REFINE )); [ $HOURS -gt 24 ] && HOURS=24

for g in $GRIDS; do
qsub <<EOF
#!/bin/bash -l
#PBS -N lt_${g}_r$REFINE
#PBS -A UMCP0023
#PBS -q casper
#PBS -l select=1:ncpus=36:mem=100GB
#PBS -l walltime=$HOURS:00:00
#PBS -o $LOGS/${g}_r$REFINE.pbs.out
#PBS -e $LOGS/${g}_r$REFINE.pbs.err

module --force purge
module load ncarenv conda
conda activate npl
cd $REPO/linear_theory
python -u stability_grid.py $g 0 1 $OUT 36 $REFINE > $LOGS/${g}_r$REFINE.log 2>&1 &&
python -u merge_grid.py $g $OUT >> $LOGS/${g}_r$REFINE.log 2>&1
EOF
done
echo "output: $OUT  (plot with: BBL_LT_DIR=$OUT python figures/<script>.py)"
