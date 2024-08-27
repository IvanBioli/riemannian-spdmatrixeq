function plot_manopt_results(results, names, fixed_rank, titletext, what_runs, what_zoom, savename, save_flag, paperflag, what_plot)
% plot_manopt_results - Plots the results of Riemannian optimization.

% Using latex
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

if nargin < 8
    save_flag = true;
end
if nargin < 9
    paperflag = false;
end
if nargin < 10
    what_plot = ones(1, 4);
end
if isstring(fixed_rank)
    rankstr = "_" + fixed_rank;
else
    rankstr = "_r" + num2str(fixed_rank);
end

savename_new = savename + rankstr;
if save_flag
    mkdir("../report/figures/" + savename_new)
end

% Plots
if what_plot(1); plotfield("cost", '$f(X_k)$', @plot); end
if what_plot(2); plotfield("gradnorm", '$||\mathrm{grad} f(X_k)||$', @semilogy); end
if what_plot(3); plotfield("res_norm_rel", 'Relative residual', @semilogy); end
if isstring(fixed_rank)
    if what_plot(4); plotfield("rank", 'Rank', @plot); end
end

% Save results
if save_flag
    save("results/manopt_" + savename_new  + ".mat", ...
        "results", "names", "fixed_rank", "titletext", "what_runs", ...
        "what_zoom", "savename")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fig = plotfield(fieldname, axisname, plotfun)
        % plotfield - Plots and saves the figure regarding one field of the info structure.
        
        % Define the color map
        if isempty(what_zoom)
            zoom_vec = [0];
        else
            zoom_vec = [0,1];
        end
        for zoom = zoom_vec
            if ~all((what_runs == 0) | (what_runs == 1))
                idx_j = what_runs;
                zoomstr = "";
            elseif zoom && all((what_runs == 0) | (what_runs == 1))
                idx_j = find(what_runs & what_zoom);
                zoomstr = "_zoom";
            else
                idx_j = find(what_runs);
                zoomstr = "";
            end
            fig = figure();
            t = tiledlayout(1, 2,'TileSpacing','Compact');
            nexttile % x-axis is iteration number
            for j = idx_j
                info = results{j};
                name = names{j};
                if isnumeric(fixed_rank)
                    name = erase(name, ", $r=" + num2str(fixed_rank) + "$");
                end
                if isfield(info, fieldname)
                    lineopts = lineopts_map(name, find(j == idx_j, 1), paperflag);
                    plotfun([info.iter], [info(:).(fieldname)], lineopts{:});
                    hold on
                end
            end
            xlabel('Iteration number $k$');
            nexttile % x-axis is time
            for j = idx_j
                info = results{j};
                name = names{j};
                if isnumeric(fixed_rank)
                    name = erase(name, ", $r=" + num2str(fixed_rank) + "$");
                end
                if isfield(info, fieldname)
                    lineopts = lineopts_map(name, find(j == idx_j, 1), paperflag);
                    plotfun([info.time], [info(:).(fieldname)], lineopts{:});
                    hold on
                end
            end
            xlabel("Time [s]")
            txt = ylabel(t, axisname); txt.Interpreter= 'latex';
            if ~isempty(titletext)
                txt = title(t, titletext); txt.Interpreter= 'latex';
            end
            set_figsize_legend(sum(idx_j, 'all'), paperflag)
            
            if save_flag
                exportgraphics(fig, "../report/figures/" + savename_new + "/manopt_" + fieldname + zoomstr + ".pdf")
            end
        end
    end

end