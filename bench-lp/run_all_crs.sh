if [ "$#" -ne 4 ]; then
  echo "Script for running Crossover for Primal-Dual Methods via COPT in batch mode..."
  echo "Usage: <dataset directory> <sol> <time limit for each instance> <crossover script>" 1>&2
  echo "Example: zsh run_all_crs_abip.sh data/netlib solution 1000 crossover_by_copt.py"
  exit -1
fi

# test set
input=$1
result=$2
timelimit=$3
# crossover
crs_script=$4

cross=$result-cross
mkdir -p $cross



# no need to run scs at 1e-8
for f in $(/bin/ls $input); do
  # echo $f
  ff=$(basename -s .mps.gz $f)
  primal=$(/bin/ls $result/${ff}*primal.txt)
  dual=$(/bin/ls $result/${ff}*dual.txt)
  cmd="timeout $timelimit python -u $crs_script $input/$f $primal $dual $cross &> $cross/$ff.log"
  echo $cmd
done
