clear;
clc;
close all;

%% TCGA input files
files = {
    'KICH_SFVA41.xlsx'
    'KIRC_SFVA41.xlsx'
    'KIRP_SFVA41.xlsx'
    'LIHC_SFVA41.xlsx'
};

labels = {
    'KICH'
    'KIRC'
    'KIRP'
    'LIHC'
};

numPathways = 41;

for i = 1:numel(files)

    fprintf('\nReading: %s\n', files{i});

    [pathways, scores] = read_FVA_xlsx(files{i});

    [scores, idx] = sort(scores, 'ascend');
    pathways = pathways(idx);

    n = min(numPathways, numel(scores));

    scores = scores(1:n);
    pathways = pathways(1:n);

    %% Narrow figure
    fig = figure( ...
        'Position', [100, 50, 520, 950], ...
        'Color', 'w');

    ax = axes( ...
        'Parent', fig, ...
        'Position', [0.57, 0.12, 0.10, 0.80]);

    imagesc(ax, scores(:));

    %% White-to-red
    nColors = 256;

    whiteRed = [ ...
        ones(nColors,1), ...
        linspace(1,0,nColors)', ...
        linspace(1,0,nColors)' ...
    ];

    colormap(ax, whiteRed);
    caxis(ax, [0.1 0.4]);

    set(ax, ...
        'YTick', 1:n, ...
        'YTickLabel', cellstr(pathways), ...
        'XTick', 1, ...
        'XTickLabel', labels{i}, ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', 8, ...
        'Box', 'on', ...
        'LineWidth', 0.5);

    title(ax, 'FVA similarity', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    draw_cell_grid(ax, n);

    cb = colorbar(ax);

    set(cb, ...
        'Position', [0.72, 0.12, 0.045, 0.80], ...
        'Limits', [0.1 0.4], ...
        'Ticks', [0.1 0.2 0.3 0.4], ...
        'FontSize', 8);

    ylabel(cb, 'score');

    outName = sprintf( ...
        'FVA_similarity_heatmap_TCGA_%s', labels{i});

    saveas(fig, [outName '.svg']);
    saveas(fig, [outName '.pdf']);

    fprintf('Created: %s.svg\n', outName);
    fprintf('Created: %s.pdf\n', outName);

    close(fig);
end

%% =========================================================
% FUNCTIONS
% =========================================================

function draw_cell_grid(ax, nRows)

    hold(ax, 'on');

    for y = 0.5:1:(nRows + 0.5)
        line(ax, [0.5 1.5], [y y], ...
            'Color', [0.45 0.45 0.45], ...
            'LineWidth', 0.25);
    end

    line(ax, [0.5 0.5], [0.5 nRows + 0.5], ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 0.4);

    line(ax, [1.5 1.5], [0.5 nRows + 0.5], ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 0.4);

    hold(ax, 'off');
end

function [pathways, scores] = read_FVA_xlsx(filename)

    if ~exist(filename, 'file')
        error('Input file not found: %s', filename);
    end

    [~, ~, raw] = xlsread(filename);

    if isempty(raw)
        error('File is empty or could not be read: %s', filename);
    end

    if size(raw,2) < 2
        error('File %s must contain at least two columns.', filename);
    end

    header = raw(1,:);
    pathwayCol = [];

    for j = 1:numel(header)

        headerName = strtrim(char(string(header{j})));

        if strcmpi(headerName, 'Pathways') || ...
           strcmpi(headerName, 'Pathway') || ...
           strcmpi(headerName, 'subSys') || ...
           strcmpi(headerName, 'Subsystem')

            pathwayCol = j;
            break;
        end
    end

    if isempty(pathwayCol)
        pathwayCol = 1;
    end

    pathwayCells = raw(2:end, pathwayCol);
    scoreCells = raw(2:end, 2);

    pathways = strings(numel(pathwayCells),1);
    scores = nan(numel(scoreCells),1);

    for r = 1:numel(pathwayCells)

        if isempty(pathwayCells{r})
            pathways(r) = "";
        else
            pathways(r) = string(pathwayCells{r});
        end

        value = scoreCells{r};

        if isnumeric(value) && isscalar(value)
            scores(r) = value;
        else
            scores(r) = str2double(string(value));
        end
    end

    valid = ...
        ~isnan(scores) & ...
        strlength(strtrim(pathways)) > 0;

    pathways = pathways(valid);
    scores = scores(valid);

    if isempty(scores)
        error('No valid scores found in %s.', filename);
    end
end