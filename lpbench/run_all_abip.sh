# /*
#  * Created Date: Monday, October 18th 2021, 11:29:23 am
#  * Author: C. Zhang
#  *
#  * Copyright (c) 2022
#  */

if [ "$#" -ne 5 ]; then
  echo "Script for generating cmds to run ABIP via matlab in batch mode..."
  echo "the commands save to ./cmd.sh"
  echo "Usage: <dataset directory> <output_directory> <time limit for each instance> <precision> <direct:0 / indirect:1>" 1>&2
  exit -1
fi

# params
wdir=$(pwd)

echo "working directory ${wdir}"

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
usedirect=$5

abipsrc='abip'
abipname=$output/${abipsrc}_${usedirect}_1e-$precision

# save to a directory
mkdir -p $abipname
if [ -f cmd.sh ]; then
  rm cmd.sh
fi

for f in $(/bin/ls $input); do
  ff=$(basename -s .mps.gz $f)

  mat_cmd="addpath matlab; addpath lpbench; test_one_abip('$input/$f', '$abipname', 1e-$precision, $timelimit); exit;"
  full_cmd="nohup timeout $timelimit \
    $MATLAB_HOME/bin/matlab -nodesktop -nodisplay -nosplash -noFigureWindows -r \
    \"${mat_cmd}\""
  echo $full_cmd
  echo $full_cmd >>cmd.sh
done
