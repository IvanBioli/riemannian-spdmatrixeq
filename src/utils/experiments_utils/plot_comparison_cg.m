function plot_comparison_cg(results, names, manifold_rank, titletext, savename, save_flag, paperflag)
% plot_comparison_cg - Plots the comparison of CG with truncation,
% fixed-rank Riemannian optimization, rank-adaptive Riemannian
% optimization, and possibly other optimization methods for multiterm
% linear matrix equations (e.g. MultiRB)

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

fig = figure();
t = tiledlayout(1, 3,'TileSpacing','Compact');

% Plot relative residual
nexttile
xaxis_name = "iter"; xaxis_label = "Iteration number $k$";
plotfield("res_norm_rel", 'Relative residual', @semilogy);
nexttile
xaxis_name = "time"; xaxis_label = "Time [s]";
plotfield("res_norm_rel", [], @semilogy);
set_figsize_legend(numel(names), paperflag)

% Plot rank
nexttile
xaxis_name = "iter"; xaxis_label = "Iteration number $k$";
plotfield("rank", 'Rank $r$', @plot);
yline(manifold_rank,'--',"$r = "+ num2str(manifold_rank) + "$", 'Interpreter','latex','HandleVisibility','off');
% Rescale fontsize
fontsize(gcf, 10, "points")
% Rescale markersize
h = findall(gcf, 'Type', 'Line');
for k = 1:length(h)
    h(k).MarkerSize = 5;
end

% Add title and save
if ~isempty(titletext)
    txt = title(t, titletext); txt.Interpreter= 'latex';
end
if save_flag
    exportgraphics(fig, "../report/figures/comparison_CGvsManopt_" + savename + ".pdf")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotfield(fieldname, axisname, plotfun)
        % plotfield - Plots and saves the figure regarding one field of the info structure.
        
        % Plots
        count = 0;
        for i = 1:length(names)
            info = results{i}; name = names{i};
            if ~isempty(info)
                count = count + 1;
                if isfield(info, fieldname)
                    name = erase(name, " (exact)");
                    lineopts = lineopts_map(name, count, paperflag);
                    plotfun([info(:).(xaxis_name)], [info(:).(fieldname)], lineopts{:});
                    hold on
                    if contains(name, "RRAM")
                        idx = find([info(1:end-1).rank] < [info(2:end).rank])+1;
                        x = [info(idx).(xaxis_name)];
                        y = [info(idx).(fieldname)];
                        scatter(x,y,10,"red","filled", 'HandleVisibility','off')
                        hold on
                    end
                end
            end
        end
        xlabel(xaxis_label)
        
        % Cut time plot if too long
        if xaxis_name == "time"
            times = [];
            for i = 1:length(names)
                info = results{i};
                if ~isempty(info)
                    times = [times, info(end).time];
                end
            end
            times = sort(times);
            if (1.5 * times(end-1)) < times(end)
                xl = xlim;
                xlim([xl(1), 1.5 * times(end-1)])
            end
            if (4 * times(end-2)) < times(end)
                xl = xlim;
                xlim([xl(1), 4 * times(end-2)])
            end
        end
        
        % Title and legend saving
        if ~isempty(axisname)
            txt = ylabel(axisname); txt.Interpreter= 'latex';
        end
    end
end