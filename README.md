<h1 align="center">
  <a href="https://github.com/IvanBioli/bioli_masters_thesis">
    <img src="report/configs/EPFL_Logo_184X53.svg" alt="EPFL Logo" width="200" height="100" hspace="20">
    <img src="report/configs/Logo_UNIPI.svg" alt="UniPi Logo" width="200" height="100" hspace="20">
  </a>
</h1>


# Preconditioned Low-Rank Riemannian Optimization for Symmetric Positive Definite Linear Matrix Equations
Code accompanying the paper "Preconditioned Low-Rank Riemannian Optimization for Symmetric Positive Definite Linear Matrix Equations".

## Reproducing Experiments

### Prerequisites

The code is written in MATLAB, R2024a. Ensure you have the appropriate MATLAB version installed.

### Installation

Clone the repository, including the Manopt submodule, using the following commands:
```
git clone https://github.com/IvanBioli/riemannian-spdmatrixeq.git
cd riemannian-spdmatrixeq
git submodule update --init 
```

### Running Experiments

To replicate all experiments detailed in the Paper, follow these steps:

1. Run the MATLAB script `src\PAPER_run_experiments.m`.

2. Allow the script to execute fully. Ignore temporary plots generated and deleted during execution; relevant plots will be generated at the end. Results are saved in `.mat` files in the folder `src\results`. All plots are saved in the folder `report\figures` and can be produced after one complete execution by running the last section of `src\PAPER_run_experiments.m.`

## Repository description
A detailed documentation can be found in each file. The repository is structured as follows.

- `report`
  - `configs` 
  - `figures`: figure folder
- `src`: code folder
  - `PAPER_run_experiments.m`: **main script** reproducing experiments and figures in the paper
  - `baseline_solvers`: non-Riemannian solvers for SPD linear matrix equations
  - `examples`: code and data for the numerical tests in the paper
    - `pdes_FD`: finite difference discretization of 2D PDEs on square domain
    - `rail_problem`: modified Bilinear Rail problem
    - `stoch_diff`: stochastic Galerkin matrix equations
  - `init.m`: adds necessary folders to MATLAB's path
  - `manifold_optimization`: code for Riemannian solvers (including manifolds and preconditioners)
  - `manopt`: fork of Manopt's repository
  - `paper_experiments`: code for setting up the numerical experimens in the paper
  - `results`: folder where results are saved in `.mat` format
  - `utils`: code utilities

## Authors
- Ivan Bioli *(Università di Pavia)*
- Daniel Kressner *(EPFL)*
- Leonardo Robol *(Università di Pisa)*
