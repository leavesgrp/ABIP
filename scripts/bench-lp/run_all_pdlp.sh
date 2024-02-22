# /*
#  * Created Date: Monday, October 18th 2021, 11:29:23 am
#  * Author: C. Zhang
#  *
#  * Copyright (c) 2022
#  */

if [ "$#" -ne 4 ]; then
  echo "Script for running Google PDLP in batch mode..."
  echo "Usage: <dataset directory> <output directory> <time limit for each instance> <eps>" 1>&2
  exit -1
fi

# params
pdhgsrc=google-pdlp

if [ ! -d $pdhgsrc ]; then
  echo "pdhg not set up here!"
  echo "link your pdgh directory to './google-pdlp' "
fi

# test set
set=$1
output=$2
timelimit=$3
precision=$4
eps=1e-$precision

name=pdhg_1e-${precision}
phdg6=$output/$name

mkdir -p $phdg6
rm cmd.pdlp.sh
for f in $(/bin/ls $set); do
    ff=$(basename -s .mps.gz $f)
    cmd="export OPENBLAS_NUM_THREADS=4; nohup julia --project=$pdhgsrc/scripts $pdhgsrc/scripts/solve_qp.jl \
      --time_sec_limit $timelimit \
      --relative_optimality_tol $eps --eps_primal_infeasible $eps --eps_dual_infeasible $eps \
      --method pdhg \
      --output_dir $phdg6 \
      --verbosity 4 \
      --redirect_stdio false \
      --instance_path $set/$f &> $phdg6/$ff.log "
    echo $cmd
    echo $cmd >> cmd.pdlp.sh
done

