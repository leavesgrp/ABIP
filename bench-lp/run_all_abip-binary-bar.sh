# /*
#  * Created Date: Monday, October 18th 2021, 11:29:23 am
#  * Author: C. Zhang
#  *
#  * Copyright (c) 2022
#  */

if [ "$#" -ne 4 ]; then
  echo "Script for generating cmds to run ABIP via binary in batch mode..."
  echo "the commands save to ./cmd.sh"
  echo "Usage: <dataset directory> <output_directory> <time limit for each instance> <eps>" 1>&2
  exit -1
fi

# params
wdir=$(pwd)

# echo "working directory ${wdir}"

# if [ -z $abipsrc ]; then
#   echo "'abipsrc' unset!"
#   echo "set your abip directory first 'abip-lp' "
#   exit
# fi

# if [ ! -d $abipsrc ]; then
#   echo "abip not set up here!"
#   echo "link your abip directory here, maybe './abip-lp' "
#   exit
# fi

# test set
input=$1
output=$2
timelimit=$3
precision=$4
finalprecision=1e-8

abipsrc=abipc-$precision-bar
abipname=$output/${abipsrc}

# save to a directory
mkdir -p $abipname
if [ -f cmd.sh ]; then
  rm cmd.sh
fi

for f in $(/bin/ls $input/*.mps); do
  ff=$(basename -s .mps $f)

  mat_cmd="bin/abip $input/$ff.mps $timelimit 100000 10000000 20 $finalprecision 1e-$precision 5 1 $abipname/$ff"
  full_cmd="nohup timeout $timelimit \
    ${mat_cmd} &> $output/${abipsrc}/$ff.log "
  echo $full_cmd
done
