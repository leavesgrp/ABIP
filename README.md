# ABIP
## ADMM-based Interior Point Method

ABIP is a new framework that applies alternating direction method of multipliers (ADMM) to implement interior point method (IPM) for solving large-scale linear programs (LP).

ABIP is developed by **[Tianyi Lin (https://github.com/tyDLin)](https://github.com/tyDLin)** and the is currently maintained by **[LEAVES](leaves.shufe.edu.cn)**  optimization software platform

## Version
Current version of ABIP is 1.00

## Reference & Performance

**Lin, Tianyi, et al. "An ADMM-based interior-point method for large-scale linear programming." Optimization Methods and Software (2020): 1-36.**

The paper is available at 
[https://arxiv.org/pdf/1805.12344.pdf](https://arxiv.org/pdf/1805.12344.pdf) with implementation details as well as reports on numerical experiments.

## Interfaces

ABIP currently supports C and MATLAB interface. 

### C
The source code of ABIP is available in `ABIP-LP\abip` directory


### MATLAB
To compile .mex files for MATLAB, run

> `make_abip.m`

in the `ABIP-LP\` directory.

