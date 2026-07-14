%% visualize_FVA_TCGA.m

clear;
clc;
close all;

inputDir = fullfile( ...
    'results','analysis','FVA','visualization_inputs');

outputDir = fullfile( ...
    'results','analysis','FVA','figures');

if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

files = {
    fullfile(inputDir,'KICH_SFVA41.xlsx')
    fullfile(inputDir,'KIRC_SFVA41.xlsx')
    fullfile(inputDir,'KIRP_SFVA41.xlsx')
    fullfile(inputDir,'LIHC_SFVA41.xlsx')
};

labels = {
    'KICH'
    'KIRC'
    'KIRP'
    'LIHC'
};

numPathways = 41;

for i = 1:numel(files)

    fprintf('\nReading: %s\n',files{i});

    [pathways,scores] = read_FVA_xlsx(files{i});

    [scores,idx] = sort(scores,'ascend');
    pathways = pathways(idx);

    n = min(numPathways,numel(scores));

    scores = scores(1:n);
    pathways = pathways(1:n);

    fig = figure( ...
        'Position',[100,50,520,950], ...
        'Color','w');

    ax = axes( ...
        'Parent',fig, ...
        'Position',[0.57,0.12,0.10,0.80]);

    imagesc(ax,scores(:));

    nColors = 256;

    whiteRed = [ ...
        ones(nColors,1), ...
        linspace(1,0,nColors)', ...
        linspace(1,0,nColors)' ...
    ];

    colormap(ax,whiteRed);
    caxis(ax,[0.1 0.4]);

    set(ax, ...
        'YTick',1:n, ...
        'YTickLabel',cellstr(pathways), ...
        'XTick',1, ...
        'XTickLabel',labels{i}, ...
        'TickLabelInterpreter','none', ...
        'FontSize',8, ...
        'Box','on', ...
        'LineWidth',0.5);

    title(ax,'FVA similarity', ...
        'FontSize',10, ...
        'FontWeight','bold');

    draw_cell_grid(ax,n);

    cb = colorbar(ax);

    set(cb, ...
        'Position',[0.72,0.12,0.045,0.80], ...
        'Limits',[0.1 0.4], ...
        'Ticks',[0.1 0.2 0.3 0.4], ...
        'FontSize',8);

    ylabel(cb,'score');

    outName = fullfile( ...
        outputDir, ...
        ['FVA_similarity_heatmap_TCGA_' labels{i}]);

    saveas(fig,[outName '.svg']);
    saveas(fig,[outName '.pdf']);

    fprintf('Created: %s.svg\n',outName);
    fprintf('Created: %s.pdf\n',outName);

    close(fig);
end

disp('TCGA FVA heatmaps finished');

%% Local functions