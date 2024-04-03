# An Enhanced ADMM-based Interior Point Method for Linear and Conic Optimization (ABIP)

This repository hosts the ongoing development of the ABIP framework, an extension of the work published in the [INFORMS Journal on Computing (IJoC)](https://pubsonline.informs.org/journal/ijoc). This project is under continuous development beyond the version documented in the paper titled "An Enhanced ADMM-based Interior Point Method for Linear and Conic Optimization" available at [https://doi.org/10.1287/ijoc.2023.0017](https://doi.org/10.1287/ijoc.2023.0017).

## Citing This Work

If this repository contributes to your research or project, please cite both the original paper and this ongoing development version using their respective DOIs.

- Paper DOI: [https://doi.org/10.1287/ijoc.2023.0017](https://doi.org/10.1287/ijoc.2023.0017)
- Repository DOI: [https://doi.org/10.1287/ijoc.2023.0017.cd](https://doi.org/10.1287/ijoc.2023.0017.cd)

Below is the BibTeX entry for citing this repository:


## Description

ADMM-based Interior Point Method for Linear and Conic Programming **(ABIP)** is a new framework that applies alternating direction method of multipliers (ADMM) to implement interior point method (IPM) for solving large-scale linear programs (LP) and conic programs (QCP).

ABIP (LP part) was initially developed by **[Tianyi Lin (https://github.com/tyDLin)](https://github.com/tyDLin)** and is currently maintained by **[LEAVES](https://github.com/leavesgrp)** optimization software platform. 

The original ABIP follows the following paper:

```
@article{lin_admm-based_2021,
	title = {An {ADMM}-based interior-point method for large-scale linear programming},
	volume = {36},
	number = {2-3},
	journal = {Optimization Methods and Software},
	author = {Lin, Tianyi and Ma, Shiqian and Ye, Yinyu and Zhang, Shuzhong},
	year = {2021},
	note = {Publisher: Taylor \& Francis},
	pages = {389--424},
}
```

### Version
Current version of ABIP is 2.0. The C interface is planned for 3.0 release.

### Installation

We provide the Matlab interface. 
To install ABIP, install [`Intel OneAPI MKL`](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit) first.

We suggest that the user set the environment variables by using the script provided by OneAPI.
It is typically located at `path-to-oneapi/setvars.sh`.

For example, on Ubuntu, the root path of OneAPI is `/opt/intel/oneapi`, then you can run the following,

```bash
cd to_the_root_directory_of_abip
source /opt/intel/oneapi/setvars.sh       
``` 

When you finish the setups, use `install.m` script in the root directory


## Basic Usage

### Basic usages for ABIP-LP

ABIP operates on  the standard form LP, i.e.,

$$
\begin{aligned}
\min ~ &c^Tx \\
\text{s.t.} ~&Ax = b\\
&x\ge 0
\end{aligned}
$$

and accepts standard `sedumi` format defined using $A, b, c, K$.

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


For convenience, we provide scripts to proceed the standard form LP, reformulates the problem to the appropriate format, see [Replicating](#replicating).

### Basic usages of ABIP-QCP

ABIP operates on the following quadratic cone programming (QCP), i.e.,

$$
\begin{aligned}
\min ~ &\frac{1}{2}x^TQx + c^Tx \\
\text{s.t.} ~&Ax = b\\
&x\in \mathcal{K}
\end{aligned}
$$

where $Q \in \mathbb{S}_{+}^n, c \in \mathbb{R}^n, b \in \mathbb{R}^m, A \in \mathbb{R}^{m \times n}$, and $\mathcal{K}$ is a closed convex cone.

To call it, use

```
[x, y, s, info] = abip(data, K, params)
```
ABIP-QCP shares most of the parameters with ABIP-LP, except for the convex cone K, currently we support following cones:
|   Parameter   | Explanation                                            |
| :-----------: | :----------------------------------------------------- |
|    K.q    | array of second-order cone constraints                                    |
|   K.rq   | array of rotated second-order cone constraints       |
|      K.f      | length of free cone     |
| K.z | length of zero cone                                |
| K.l  | length of LP cone                                  |

**Note**: columns of data matrix A must be specified in the order above.

## Results

The result files can be found in `results`

## Replicating

To replicate the results in the paper, please refer to the README [here](scripts/README.md)

## Ongoing Development

This code is being developed on an on-going basis at the author's
[Github site](https://github.com/leavesgrp/ABIP).

## Support

For support in using this software, submit an
[issue](https://github.com/leavesgrp/ABIP/issues/new).

## License

### Project License

This project is licensed under the MIT License, as detailed in the `LICENSE` file.

### External Dependencies and Licenses

All external dependencies are consolidated in `src/external/`, with each dependency's license:

- **AMD (SuiteSparse)**: BSD license allows broad usage with minimal restrictions. The license file is included in its respective directory.
- **LDL & CSparse (SuiteSparse)**: Under GNU LGPL, permitting use and modification provided changes are shared under LGPL. License files are placed within their specific directories.
- **QDLDL (OSQP)**: Apache License 2.0 supports use and modification, requiring clear documentation of changes and the original copyright notice. The license file is located in the relevant directory.

### Contributions

Contributors agree that their contributions will be licensed under the MIT License. We encourage contributors to review the licenses of external dependencies housed in `src/external/` to ensure compliance with their terms.


---

This README section is provided for informational purposes only and does not constitute legal advice. If you have any questions or concerns about the licensing terms or compliance, especially regarding third-party dependencies, please consult with a legal expert.
