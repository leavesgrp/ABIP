## ABIP
### ADMM-based Interior Point Method for Linear and Conic Programming

ABIP is a new framework that applies alternating direction method of multipliers (ADMM) to implement interior point method (IPM) for solving large-scale linear programs (LP) and conic programs (QCP.

ABIP(LP) was initially developed by **[Tianyi Lin (https://github.com/tyDLin)](https://github.com/tyDLin)** and is currently maintained by **[LEAVES](https://leaves.shufe.edu.cn)**  optimization software platform. 

### Version
Current version of ABIP is 2.0

v2.0 (ABIP+)

1. Added support for QCP
2. Improved LP solver

v1.0

1. ABIP for large scale linear programming problems

### Reference

- Lin, Tianyi, et al. "An ADMM-based interior-point method for large-scale linear programming." Optimization Methods and Software (2020): 1-36.
- Deng, Qi, et al. "New Developments of ADMM-based Interior Point Methods for Linear Programming and Conic Programming." *arXiv preprint arXiv:2209.01793* (2022).

### Installation

To install ABIP, use `install.m` script in the main direction.

For QCP solver,  user needs to specify MKL path in 

```
mkl_path = ''; % '/opt/intel/oneapi/mkl'
```

### How to call ABIP

ABIP accepts standard sedumi format defined using $A, b, c, K$.

To call it, use

```
[x, y, s, info] = abip(data, K, params)
```

|   Parameter   | Explanation                                            |
| :-----------: | :----------------------------------------------------- |
|    verbose    | If log is turned on                                    |
|   normalize   | Whether to perform data scaling before the solve       |
|      pcg      | Whether to use iterative solver for linear systems     |
| max_admm_iter | Maximum ADMM iteration                                 |
| max_ipm_iter  | Maximum IPM iteration                                  |
|   timelimit   | Time limit                                             |
|      tol      | Relative tolerance of convergence                      |
|    solver     | Choose which solver to use, -1: automatic 1: force QCP |

### Current developers

**ABIP-LP**

Chenyu Xue (Phd. student at RIIS, SUFE. E-mail: xcy2721d@gmail.com)

**ABIP-QCP**

Jinsong Liu (Phd. student at SUFE. Email: liujinsong@163.shufe.edu.cn)

Tianhao Liu (Phd. student at RIIS, SUFE. Email: liu.tianhao@163.sufe.edu.cn)
