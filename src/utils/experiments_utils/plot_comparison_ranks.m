function plot_comparison_ranks(results, names, ranks, titletext, savename, save_flag, paperflag, field)
% plot_comparison_ranks - Plots the comparison of fixed-rank Riemannian
% optimization with different (fixed) ranks.

% Using latex
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

if nargin < 6
    save_flag = true;
end
if nargin < 7
    paperflag = false;
end
if nargin < 8
    field.name = "gradnorm";
    field.label = "Gradient norm";
end

fig = figure();
t = tiledlayout(1, 2,'TileSpacing','Compact');
xaxis_name = "iter"; xaxis_label = "Iteration number $k$";

% Plot results first rank
nexttile
plotfield(results{1}, names{1}, ranks(1), field.name, field.label, @semilogy);
txt = title("$ r = "+ num2str(ranks(1)) + "$"); txt.Interpreter= 'latex';
% Plot results second rank
nexttile
plotfield(results{2}, names{2}, ranks(2), field.name, field.label, @semilogy);
txt = title("$ r = "+ num2str(ranks(2)) + "$"); txt.Interpreter= 'latex';
set_figsize_legend(numel(names), paperflag)
% Add legend and title
if ~isempty(titletext)
    txt = title(t, titletext); txt.Interpreter= 'latex';
end
% Save
if save_flag
    exportgraphics(fig, "../report/figures/comparison_r" + num2str(ranks(1)) + ...
        "_r" + num2str(ranks(2)) + "_" + savename + ".pdf")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotfield(results, names, rank, fieldname, axisname, plotfun)
        % plotfield - Plots and saves the figure regarding one field of the info structure.
        
        % Plots
        for i = 1:length(names)
            info = results{i};
            if ~isempty(info)
                name = erase(names{i}, ", $r=" + num2str(rank) + "$");
                if isfield(info, fieldname)
                    lineopts = lineopts_map(name, i, paperflag);
                    plotfun([info(:).(xaxis_name)], [info(:).(fieldname)], lineopts{:});
                    hold on
                end
            end
        end
        xlabel(xaxis_label)
        
        % Title and legend saving
        txt = ylabel(axisname); txt.Interpreter= 'latex';
    end
end