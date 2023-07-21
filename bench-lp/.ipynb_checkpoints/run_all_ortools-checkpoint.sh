# /*
#  * Created Date: Monday, October 18th 2021, 11:29:23 am
#  * Author: C. Zhang
#  *
#  * Copyright (c) 2022
#  */

if [ "$#" -ne 4 ]; then
  echo "Script for running Google PDLP in batch mode..."
  echo "Usage: <dataset directory> <output> <time limit for each instance> <eps>" 1>&2
  exit -1
fi

wdir=$(pwd)
# test set
data=$1
output=$2
timelimit=$3
precision=$4
eps=1e-$precision

name=pdhg_cpp_1e-${precision}
phdg=$output/$name

mkdir -p $phdg

cd $wdir/$pdhgsrc
for f in $(/bin/ls $data/*.mps); do

  ff=$(basename -s .mps $f)
  cmd="python abip-lp/bench-lp/pdlp_solve.py \
    --file $data/$ff.mps \
    --tol  $eps \
    --output $phdg
    &> $phdg/$ff.log"

  echo $cmd 

done
