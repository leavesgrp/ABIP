# /*
#  * Created Date: Monday, October 18th 2021, 11:29:23 am
#  * Author: C. Zhang
#  *
#  * Copyright (c) 2022
#  */

if [ "$#" -ne 2 ]; then
  echo "Script for generating cmds to run ABIP via matlab in batch mode..."
  echo "the commands save to ./cmd.sh"
  echo "Usage: <dataset directory> <output_directory> <time limit for each instance>" 1>&2
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

abipsrc='abip'
abipname=$output/

# save to a directory
mkdir -p $abipname
if [ -f cmd.sh ]; then
  rm cmd.sh
fi

for f in $(/bin/ls $input); do
  ff=$(basename -s .mps.gz $f)

  mat_cmd="addpath '/home/chuwen/workspace/abip-lp/src/abip-lp/interface'; \
  addpath '/home/chuwen/gurobi1001/linux64/matlab/'; \
  addpath matlab; addpath 'bench-lp'; save_abip_mps('$input/$f', '$abipname'); exit;"
  full_cmd="nohup timeout 3000 \
    $MATLAB_HOME/bin/matlab -nodesktop -nodisplay -nosplash -noFigureWindows -r \
    \"${mat_cmd}\" &> $output/$f.log "
  echo $full_cmd
  echo $full_cmd >>cmd.presolve.sh
done
