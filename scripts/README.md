# A Guide for Benchmarking
We assume the access to standard libraries such as the MIPLIB, Netlib, cblib and so on.

## Data Preparation
See [Preparing Data](../data/README.md).

## Benchmarking ABIP-LP
The scripts are collected at `scripts/bench-lp`. 

ABIP uses the standard form LP. For convenience we provide the following scripts:

- preprocess the data in the MPS and produce a standard form LP: `scripts/bench-lp/preprocess.m`
- saving standard form LP from a mps file `scripts/bench-lp/save_abip_mps.m`
- running ABIP with default parameters`scripts/bench-lp/test_one_abip.m`
    - the script summarizes the solutions and statistics into a json file.

Built upon this, we include other competing methods.
For PDLP install `ortools`:

```
pip install ortools==9.3.10497
```

- PDLP: `scripts/bench-lp/pdlp_solve.py`
- COPT: `scripts/bench-lp/copt_solve_lp.py`

Similarly, both outputs statistics into a json file. 

For creating the batch tests,
- use the bash scripts `run_*` that outputs the demand needed for each instance. 
- then analyze the json files using `analyze_*.py`. 

Then the results in the paper are summarized in the excel files in `results`.

### Note for pagerank
The matrices used to create the pagerank instances are acquired from https://sparse.tamu.edu/.

## Benchmarking ABIP-QCP

After installation, if you want to reproduce the result of **QCP**(**Lasso**, **SVM** and **cblib**) in the paper, you need to specify the data path in `test_cblib.m`, `test_svm.m` and `cblib_test_cosmo.py`, then run `test_all.m` and `cblib_test_cosmo.py`.

