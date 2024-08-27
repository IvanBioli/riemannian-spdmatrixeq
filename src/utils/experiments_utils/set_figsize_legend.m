function set_figsize_legend(num_elements, paperflag, singlefig)
% set_figsize_legend - Function for setting the size of the figure and the 
% position of the legend

    if nargin < 3
        singlefig = false;
    end
    if ~paperflag
        if num_elements <= 8
            lgnd = legend('NumColumns', 2);
            lgnd.Layout.Tile = 'south';
            pos = get(gcf, 'Position'); pos(3) = 1.5 * pos(3); pos(4) = 1.1 * pos(4); set(gcf, 'Position', pos)
        else
            lgnd = legend('NumColumns', 1);
            lgnd.Layout.Tile = 'east';
            pos = get(gcf, 'Position'); pos(3) = 3 * pos(3); set(gcf, 'Position', pos)
        end
    else
        % Set legent properties
        lgnd = legend('NumColumns', 2 - singlefig);
        if ~singlefig
            lgnd.Layout.Tile = 'south';
        else
            lgnd.Location = 'southoutside';
        end
        % Rescale figure
        set(gcf, 'Units', 'centimeters'); pos = get(gcf, 'Position'); 
        if ~singlefig
            pos(3) = 1.4 * pos(3); set(gcf, 'Position', pos)
            textwidth = 14.68 * 1.6; 
            pos = get(gcf, 'Position'); pos(4) = pos(4) * textwidth / pos(3); pos(3) = textwidth; set(gcf, 'Position', pos); 
        else
            pos(4) = 1.2 * pos(4); set(gcf, 'Position', pos)
            textwidth = 14.68 * 1.5 / 2.1; 
            pos = get(gcf, 'Position'); pos(4) = pos(4) * textwidth / pos(3); pos(3) = textwidth; set(gcf, 'Position', pos); 
        end
        % Rescale fontsize
        fontsize(gcf, 10, "points")
        % Rescale markersize
        h = findall(gcf, 'Type', 'Line');
        for i = 1:length(h)
            h(i).MarkerSize = 5;
        end
    end
end

