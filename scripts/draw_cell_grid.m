function draw_cell_grid(ax,nRows)
%DRAW_CELL_GRID Draw borders around heatmap cells.
%
% MATLAB R2019b compatible.

hold(ax,'on');

for y = 0.5:1:(nRows + 0.5)

    line(ax,[0.5 1.5],[y y], ...
        'Color',[0.45 0.45 0.45], ...
        'LineWidth',0.25);
end

line(ax,[0.5 0.5],[0.5 nRows + 0.5], ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.4);

line(ax,[1.5 1.5],[0.5 nRows + 0.5], ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',0.4);

hold(ax,'off');

end