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

2. Allow the script to execute fully. Ignore temporary plots generated and deleted during execution; relevant plots will be generated at the end. Results are saved in `.mat` files in the folder `src\results`. All plots are saved in the folder `report\figures` and can be produced after one complete execution by running the last section of `src\run_experiments.m.`

## Authors
- Ivan Bioli *(Università di Pavia)*
- Daniel Kressner *(École Polytechnique Fédérale de Lausanne - EPFL)*
- Leonardo Robol *(Università di Pisa)*
