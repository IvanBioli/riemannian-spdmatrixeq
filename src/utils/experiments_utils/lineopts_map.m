function lineopts = lineopts_map(name, j, paper_colors)
% lineopts_map - Defines the line options for each plot in the report, so
% that colors and linestyles are consistend throughout numerical
% experiments.

if nargin < 3
    paper_colors = false;
end

palette =[...
    0         0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    1.0000    0.2700    0.2270;...
    0.3960    0.5090    0.9920;...
    1.0000    0.8390    0.0390;...
    0    0.6390    0.6390;...
    0.7960    0.5170    0.3640];

%%%%%%%%%%%%%%%%%%%%%%%%%% USUAL MATLAB COLORMAP %%%%%%%%%%%%%%%%%%%%%%%%%%
if ~paper_colors
    if mod(j, 12) == 0
        idx_col = 12;
    else
        idx_col = mod(j,12);
    end
    col = palette(idx_col, :);
    marker = ".";
    linestyle = "-";
    lineopts = {"Color", col, "Marker", marker, "MarkerSize", 5, "DisplayName", name, "Linestyle", linestyle};
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% COLORMAP FOR REPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Color based on the optimization method
if contains(name, "CG trunc")
    if contains(name, "cap")
        col_idx = 5;
    else
        col_idx = "black";
    end
elseif contains(name, "CG (full rank)")
    col_idx = "black";
elseif contains(name, "Prec. Rich.") 
    if contains(name, "R-NLCG")
        col_idx = 10;
    elseif contains(name, "R-GD")
        col_idx = "magenta";
    end
elseif contains(name, "R-GD")
    col_idx = [.7 .7 .7];
elseif contains(name, "R-NLCG")
    if contains(name, "RRAM")
        col_idx = 3;
    else
        col_idx = 1;
    end
elseif contains(name, "R-TR Riemann. Hess")
    col_idx = 6;
elseif contains(name, "R-TR Proj. EHess")
    if contains(name, "RRAM")
        col_idx = 4;
    else
        col_idx = 2;
    end
elseif contains(name, "MultiRB")
    col_idx = 12;
else
    error("Unknown optimization method for legend")
end
if isnumeric(col_idx) && isscalar(col_idx)
    col = palette(col_idx, :);
else
    col = col_idx;
end

% Assign marker based on the preconditioner
if contains(name, "R-NLCG") || contains(name, "R-TR Proj. EHess")
    %----------------- EXAMPLE 1 -----------------%
    if contains(name, "tangADI") && (contains(name, "AXD+DXA") || contains(name, "AXD+EXB"))
        marker = "o";
    elseif (contains(name, "AX+XA") || contains(name, "AX+XB"))
        marker = "*";
    elseif (contains(name, "AXD+DXA") || contains(name, "AXD+EXB"))
        marker = ".";
        %----------------- EXAMPLE 2 -----------------%
    elseif contains(name, "K_0X") || contains(name, "K_0XG")
        marker = ".";
        %----------------- EXAMPLE 3 -----------------%
    elseif contains(name, "tangADI") && contains(name, "AXM+MXA")
        marker = "o";
    elseif ~contains(name, "tangADI") && contains(name, "AXM+MXA")
        marker = ".";
        %---------------------------------------------%
    else
        error("Unrecognized preconditioner")
    end
else
    marker = ".";
end

% Define linestyle
if contains(name, "CG") && contains(name, "full rank")
    linestyle = "--";
else
    linestyle = "-";
end

lineopts = {"Color", col, "Marker", marker, "MarkerSize", 5, "DisplayName", name, "Linestyle", linestyle};
end

