% Description: Runs all numerical experiments and produces the plots in the report
init;
clear all; close all; clc;
seed = 0; % Set seed for reproducibility
%% PDEs with separable coefficients and finite differences
PAPER_experiments_PDEs_FD;
clearvars -except seed; close all; clc;

%% Stochastic Galerkin matrix equation
PAPER_experiments_stoch_diff;
clearvars -except seed; close all; clc;

%% Rail problem
PAPER_experiments_rail;
clear all; close all; clc;

%% Produce all the plots
PAPER_plot_PDEs_FD;
PAPER_plot_stoch_diff;
PAPER_plot_rail;