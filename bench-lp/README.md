# Benchmarking ABIP

## Install

See `../install.m` for details

Ensure `ortools`:
```
pip install ortools==9.3.10497
```
## How to use 

We have a set of bash scipts:
```bash
bash run_*.sh
```

Once you run this script, the commands will output to the stdout. Then you can run the command one-by-one, or use `xargs` to run in parallel.

For convenience, we provide example outputs from these test scripts, see `bench-lp/sbin`.


## Number of Threads

For benchmarking, we keep 4 threads.

```
export MKL_NUM_THREADS=4; export NUMEXPR_NUM_THREADS=4; export OMP_NUM_THREADS=4; ...
```