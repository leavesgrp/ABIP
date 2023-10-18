# test set
input=$1
set=$2
simplex=$set/copt_sol_simplex
dual=$set/copt_sol_dual
barrier=$set/copt_sol_barrier

# clean
#mkdir -pv $simplex && mkdir $simplex
#mkdir -pv $dual && mkdir $dual
mkdir -pv $barrier

# gurobi
# python copt_solve_lp.py $set/copt_presolved $simplex ./copt.simplex.par &>$set/copt.simplex.log &
# python copt_solve_lp.py $set/copt_presolved $dual ./copt.dual.par &>$set/copt.dual.log &
python bench-lp/copt_solve_lp.py $input $barrier bench-lp/copt.barrier.par &> $set/copt.barrier.log &
