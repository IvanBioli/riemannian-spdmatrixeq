function plot_cg_results(results, names, titletext, savename, save_flag, paperflag)
% plot_cg_results - Plots the results of CG with truncation

% Using latex
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

if nargin < 5
    save_flag = true;
end
if nargin < 6
    paperflag = false;
end

plotfield("rank", 'Rank', @plot);
plotfield("res_norm_rel", 'Relative residual', @semilogy);

if save_flag
    save("results/cg_" + savename  + ".mat", "results", "names", ...
        "titletext", "savename")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fig = plotfield(fieldname, axisname, plotfun)
        % plotfield - lots and saves the figure regarding one field of the info structure.
        
        % Define colors and linestyles
        linestyles = ["-", "--", ":", "-."];
        
        % Plots
        fig = figure();
        t = tiledlayout(1, 2,'TileSpacing','Compact');
        
        nexttile % x-axis is iter
        for i = 1:size(names, 1)
            for j = 1:size(names, 2)
                info = results{i, j}; name = names{i, j};
                if ~isempty(info)
                    lineopts = lineopts_map(name, j);
                    plotfun([info.iter], [info(:).(fieldname)], lineopts{:});
                    hold on
                end
            end
        end
        xlabel("Iteration number")
        
        nexttile % x-axis is time
        for i = 1:size(names, 1)
            for j = 1:size(names, 2)
                info = results{i, j}; name = names{i, j};
                if ~isempty(info)
                    lineopts = lineopts_map(name, j);
                    plotfun([info.time], [info(:).(fieldname)], lineopts{:});
                    hold on
                end
            end
        end
        xlabel("Time [s]")
        
        % Title, legend and saving
        txt = ylabel(t, axisname); txt.Interpreter= 'latex';
        txt = title(t, titletext); txt.Interpreter= 'latex';
        set_figsize_legend(numel(names), paperflag)
        
        if save_flag
            exportgraphics(fig, "../report/figures/cg_" + savename + "_" + fieldname + ".pdf")
        end
    end
end