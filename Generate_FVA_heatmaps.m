clear;
clc;
close all;

%% =========================================================
% INPUT FILES
% ==========================================================

% RCC / 769-P cell-line files
files_RCC = {
    '7sc1.xlsx'
    '7sc2.xlsx'
    '7sc12.xlsx'
};

labels_RCC = {
    'RCC sc1 vs ctrl'
    'RCC sc2 vs ctrl'
    'RCC sc12 vs ctrl'
};

outputNames_RCC = {
    'sc1'
    'sc2'
    'sc12'
};

% HCC / Huh7 cell-line files
files_HCC = {
    'Hsc1.xlsx'
    'Hsc2.xlsx'
    'Hsc12.xlsx'
};

labels_HCC = {
    'HCC sc1 vs ctrl'
    'HCC sc2 vs ctrl'
    'HCC sc12 vs ctrl'
};

outputNames_HCC = {
    'sc1'
    'sc2'
    'sc12'
};

% TCGA files
files_TCGA = {
    'KICH_SFVA41.xlsx'
    'KIRC_SFVA41.xlsx'
    'KIRP_SFVA41.xlsx'
    'LIHC_SFVA41.xlsx'
};

labels_TCGA = {
    'TCGA KICH ON vs OFF'
    'TCGA KIRC ON vs OFF'
    'TCGA KIRP ON vs OFF'
    'TCGA LIHC ON vs OFF'
};

outputNames_TCGA = {
    'KICH'
    'KIRC'
    'KIRP'
    'LIHC'
};

%% =========================================================
% CREATE HEATMAPS
% ==========================================================

make_FVA_heatmaps( ...
    files_RCC, labels_RCC, outputNames_RCC, 'RCC_7', 40);

make_FVA_heatmaps( ...
    files_HCC, labels_HCC, outputNames_HCC, 'HCC_H', 40);

make_FVA_heatmaps( ...
    files_TCGA, labels_TCGA, outputNames_TCGA, 'TCGA', 41);

%% =========================================================
% LOCAL FUNCTION
% ==========================================================

function make_FVA_heatmaps( ...
    files, labels, outputNames, outPrefix, numPathways)

    if numel(files) ~= numel(labels) || ...
       numel(files) ~= numel(outputNames)

        error(['files, labels, and outputNames must contain ' ...
               'the same number of entries.']);
    end

    for i = 1:numel(files)

        %% Check input file
        if ~isfile(files{i})
            warning('File not found: %s', files{i});
            continue;
        end

        %% Read Excel table
        T = readtable( ...
            files{i}, ...
            'VariableNamingRule', 'preserve');

        if width(T) < 2
            warning(['File %s must contain at least two columns: ' ...
                     'pathway and similarity score.'], files{i});
            continue;
        end

        variableNames = T.Properties.VariableNames;

        %% Detect pathway column
        pathwayIndex = find( ...
            strcmpi(variableNames, 'Pathways'), 1);

        if isempty(pathwayIndex)
            pathwayIndex = find( ...
                strcmpi(variableNames, 'subSys'), 1);
        end

        if isempty(pathwayIndex)
            pathwayIndex = 1;
        end

        pathwayNames = string(T{:, pathwayIndex});

        %% Similarity scores are expected in column 2
        scores = T{:, 2};

        if ~isnumeric(scores)
            scores = str2double(string(scores));
        end

        %% Remove empty or invalid rows
        validRows = ...
            ~isnan(scores) & ...
            ~ismissing(pathwayNames) & ...
            strlength(strtrim(pathwayNames)) > 0;

        pathwayNames = pathwayNames(validRows);
        scores = scores(validRows);

        if isempty(scores)
            warning('No valid scores found in %s.', files{i});
            continue;
        end

        %% Sort ascending
        % Lowest similarity means greatest difference.
        [scoresSorted, order] = sort(scores, 'ascend');

        n = min(numPathways, numel(scoresSorted));

        topScores = scoresSorted(1:n);
        topPathways = pathwayNames(order(1:n));

        %% Create heatmap
        fig = figure( ...
            'Position', [100, 100, 700, 1000], ...
            'Color', 'w');

        imagesc(topScores(:));

        colormap(flipud(hot(256)));

        cb = colorbar;
        cb.Label.String = 'FVA similarity score';

        clim([0 1]);

        set(gca, ...
            'YTick', 1:n, ...
            'YTickLabel', topPathways, ...
            'TickLabelInterpreter', 'none', ...
            'XTick', 1, ...
            'XTickLabel', labels{i}, ...
            'XTickLabelRotation', 45, ...
            'FontSize', 10, ...
            'FontWeight', 'bold');

        title( ...
            sprintf('Top %d Different Pathways: %s', ...
                    n, labels{i}), ...
            'Interpreter', 'none');

        xlabel('Model comparison');
        ylabel('Pathways');

        %% Output filename
        outName = sprintf( ...
            'FVA_similarity_heatmap_%s_%s', ...
            outPrefix, outputNames{i});

        %% Save files
        exportgraphics( ...
            fig, [outName '.svg'], ...
            'ContentType', 'vector');

        exportgraphics( ...
            fig, [outName '.pdf'], ...
            'ContentType', 'vector');

        fprintf('\nInput:  %s\n', files{i});
        fprintf('Output: %s.svg\n', outName);
        fprintf('Output: %s.pdf\n', outName);

        close(fig);
    end
end