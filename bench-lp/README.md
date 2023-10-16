# Benchmarking ABIP

## Install

See `../install.m` for details

Ensure `ortools`:
```
pip install ortools==9.3.10497
```
## How to use 

Invoke the following command in bash:
```bash
bash run_all_abip.sh
```
to see the helper information.

An example run could be: (suppose working directory is at `ABIP/`)
```bash
bash lpbench/run_all_abip.sh \
    dir_in \
    dir_out \
    1000 \
    4 \
    0
```
where 
```
- dir_in  <a directory containing mps files>
- dir_out <the output directory to save results>
- 1000    <the timelimit of each instance>
- 4       <the precision, 4 means to solve a problem with tolerance of 0.0001>
- 0       <0-use direct method, 1-use PCG to solve linear systems>
```

Once you run this script, you will see a file `cmd.sh` in this directory. Then you can run the command one-by-one, or use `xargs` to run in parallel.


## Number of Threads

For benchmarking, we keep 4 threads,

```
export MKL_NUM_THREADS=4; export NUMEXPR_NUM_THREADS=4; export OMP_NUM_THREADS=4; ...
```