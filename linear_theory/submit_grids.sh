#!/bin/bash
# Submit the linear-theory parameter grids (corrected Eq 2.27 system) as Casper CPU jobs.
# Usage (from the repo root): [REFINE=n] [DM=step] [HOURS=h] bash linear_theory/submit_grids.sh [sgamma] [freq]
#   REFINE=n divides the grid spacing by n (default 1, the preprint grids); cost ~ n^2.
#   DM is the m-scan step (default 0.1, the preprint value); DM=0.01 costs ~10x.
# Output: $OUT/<grid>.nc, with OUT = output/linear_theory[_r<n>][_dm<step>];
# live progress in logs/linear_theory/<grid><suffix>.log

REFINE=${REFINE:-1}
DM=${DM:-0.1}
SUFFIX=""
[ "$REFINE" != 1 ] && SUFFIX="${SUFFIX}_r$REFINE"
[ "$DM" != 0.1 ] && SUFFIX="${SUFFIX}_dm$DM"
OUT=/glade/derecho/scratch/$USER/BottomBoundaryLayer/output/linear_theory$SUFFIX
REPO=$(cd "$(dirname "$0")/.." && pwd)
LOGS=$REPO/logs/linear_theory
mkdir -p "$OUT" "$LOGS"
GRIDS=${@:-sgamma freq}
HOURS=${HOURS:-8}    # walltime per grid job

for g in $GRIDS; do
qsub <<EOF
#!/bin/bash -l
#PBS -N lt_${g}$SUFFIX
#PBS -A UMCP0023
#PBS -q casper
#PBS -l select=1:ncpus=36:mem=100GB
#PBS -l walltime=$HOURS:00:00
#PBS -o $LOGS/${g}$SUFFIX.pbs.out
#PBS -e $LOGS/${g}$SUFFIX.pbs.err

module --force purge
module load ncarenv conda
conda activate npl
cd $REPO/linear_theory
python -u stability_grid.py $g 0 1 $OUT 36 $REFINE $DM > $LOGS/${g}$SUFFIX.log 2>&1 &&
python -u merge_grid.py $g $OUT >> $LOGS/${g}$SUFFIX.log 2>&1
EOF
done
echo "output: $OUT  (plot with: BBL_LT_DIR=$OUT python figures/<script>.py)"
