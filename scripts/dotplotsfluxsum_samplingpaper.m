clear;
clc;
close all;

%% =========================================================
% CLEAN MATLAB R2019b-COMPATIBLE SCRIPT
% Save this file as:
% dotplotsfluxsum_samplingpaper.m
% ==========================================================

%% =========================================================
% INPUT AND OUTPUT FOLDERS
% ==========================================================

inputFolder = [ ...
    '\\atlas.uni.lux\fstc_sysbio\0- UserFolders\' ...
    'Evelyn.GONZALEZ\Chiarareview\' ...
    'akr1a1_stats_sampling_medium'];

outputFolder = fullfile( ...
    inputFolder, ...
    'results', ...
    'paper_dotplots_matlab');

if ~exist(inputFolder, 'dir')
    error('Input folder not found:\n%s', inputFolder);
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

fprintf('Input folder:\n%s\n\n', inputFolder);
fprintf('Output folder:\n%s\n\n', outputFolder);

%% =========================================================
% INPUT FILES
% ==========================================================

fluxSumFiles = {
    'LIHC_fluxSumStatst_all_medium_4_7_VS_8.xlsx'
    'KIRC_fluxSumStatst_all_medium_4_1_VS_2.xlsx'
    'Huh7_fluxSumStatst_all_medium_4_5_VS_71500.xlsx'
    '769-P_fluxSumStatst_all_medium_4_1_VS_31500.xlsx'
};

samplingFiles = {
    'LIHC-ONOFF_FluxSamplingStatst_all_medium_5onoff7_vs_8.xlsx'
    'KIRC-ONOFF_FluxSamplingStatst_all_medium_5onoff1_vs_2.xlsx'
    'Huh7_FluxSamplingStatst_all_medium_55_vs_7_1500.xlsx'
    '769-P_FluxSamplingStatst_all_medium_51_vs_3_1500.xlsx'
};

conditions = {
    'LIHC'
    'KIRC'
    'Huh7'
    '769-P'
};

%% =========================================================
% FIGURE d: GLYCOLYSIS METABOLITES
% ==========================================================

glycolysisIDs = {
    'glc_D[c]'
    'g6p[c]'
    'f6p[c]'
    'fdp[c]'
    'g3p[c]'
    'dhap[c]'
    '13dpg[c]'
    '3pg[c]'
    '2pg[c]'
    'pep[c]'
    'pyr[c]'
    'lac_L[e]'
};

glycolysisLabels = {
    'glc_D'
    'g6p'
    'f6p'
    'fdp'
    'g3p'
    'dhap'
    '1,3dpg'
    '3pg'
    '2pg'
    'pep'
    'pyr'
    'lac_L'
};

glycolysisNames = {
    'D-Glucose'
    'Glucose 6-phosphate'
    'Fructose 6-phosphate'
    'Fructose 1,6-bisphosphate'
    'Glyceraldehyde 3-phosphate'
    'Dihydroxyacetone phosphate'
    '1,3-Bisphosphoglycerate'
    '3-Phosphoglycerate'
    '2-Phosphoglycerate'
    'Phosphoenolpyruvate'
    'Pyruvate'
    'L-Lactate'
};

glycolysisFC = readSelectedMetabolites( ...
    inputFolder, ...
    fluxSumFiles, ...
    glycolysisIDs);

makePaperDotPlot( ...
    glycolysisFC, ...
    glycolysisLabels, ...
    conditions, ...
    'Glycolysis Pathway', ...
    fullfile(outputFolder, 'Figure_d_Glycolysis_Pathway'));

writeSelectedMatrix( ...
    fullfile(outputFolder, 'Figure_d_Glycolysis_data.xlsx'), ...
    glycolysisIDs, ...
    glycolysisNames, ...
    conditions, ...
    glycolysisFC);

%% =========================================================
% FIGURE e: PYRUVATE METABOLITES
% ==========================================================

pyruvateMetIDs = {
    'mthgxl[c]'
    '12ppd_S[c]'
    '12ppd_R[c]'
    'lald_L[c]'
    'lald_L[m]'
    'lald_D[c]'
    'ac[m]'
    'aac[c]'
    'acetol[c]'
    'lgt_S[c]'
    'lgt_S[m]'
    'lac_D[c]'
    'lac_D[m]'
    'lac_L[m]'
    'pyr[m]'
    'oaa[m]'
    'acacoa[c]'
    'acacoa[m]'
    'gthrd[c]'
    'gthrd[m]'
};

pyruvateMetLabels = {
    'mthgxl (c)'
    '12ppd_S (c)'
    '12ppd_R (c)'
    'lald_L (c)'
    'lald_L (m)'
    'lald_D (c)'
    'ac (m)'
    'aac (c)'
    'acetol (c)'
    'lgt_S (c)'
    'lgt_S (m)'
    'lac_D (c)'
    'lac_D (m)'
    'lac_L (m)'
    'pyr (m)'
    'oaa (m)'
    'acacoa (c)'
    'acacoa (m)'
    'gthrd (c)'
    'gthrd (m)'
};

pyruvateFC = readSelectedMetabolites( ...
    inputFolder, ...
    fluxSumFiles, ...
    pyruvateMetIDs);

makePaperDotPlot( ...
    pyruvateFC, ...
    pyruvateMetLabels, ...
    conditions, ...
    'Pyruvate metabolism Pathway', ...
    fullfile(outputFolder, 'Figure_e_Pyruvate_Metabolites'));

writeSelectedMatrix( ...
    fullfile(outputFolder, 'Figure_e_Pyruvate_metabolites_data.xlsx'), ...
    pyruvateMetIDs, ...
    pyruvateMetLabels, ...
    conditions, ...
    pyruvateFC);

%% =========================================================
% FIGURE f: PYRUVATE REACTIONS
% ==========================================================

selectedReactions = {
    'MGSA'
    'MGSA2'
    'ALR2'
    'ALR3'
    'ALCD21_L'
    'ALCD22_L'
    'LCADI'
    'GLYOX'
    'LCADI_D'
    'HMR_8501'
    'ALCD22_D'
    'ALCD21_D'
    'LCADIm'
    'LGTHL'
    'PPDOy'
    'ME2'
    'LDH_Lm'
    'HMR_3855'
};

reactionFC = readSelectedReactions( ...
    inputFolder, ...
    samplingFiles, ...
    selectedReactions, ...
    'Pyruvate metabolism');

makePaperDotPlot( ...
    reactionFC, ...
    selectedReactions, ...
    conditions, ...
    'Pyruvate metabolism Pathway (reactions)', ...
    fullfile(outputFolder, 'Figure_f_Pyruvate_Reactions'));

writeSelectedMatrix( ...
    fullfile(outputFolder, 'Figure_f_Pyruvate_reactions_data.xlsx'), ...
    selectedReactions, ...
    selectedReactions, ...
    conditions, ...
    reactionFC);

fprintf('\nFinished.\n');
fprintf('Files saved in:\n%s\n', outputFolder);

%% =========================================================
% LOCAL FUNCTIONS
% IMPORTANT: keep every function below all script commands.
% ==========================================================

function result = readSelectedMetabolites(folder, files, selectedIDs)

    nFiles = numel(files);
    nMetabolites = numel(selectedIDs);

    result = nan(nFiles, nMetabolites);

    for f = 1:nFiles

        filename = fullfile(folder, files{f});

        if ~exist(filename, 'file')
            error('File not found:\n%s', filename);
        end

        fprintf('Reading flux-sum file:\n%s\n', filename);

        [headers, raw] = readExcelRaw(filename);

        idCol = findColumn(headers, { ...
            'Row', ...
            'mets', ...
            'Metabolite_ID', ...
            'MetaboliteID', ...
            'metID', ...
            '...1'});

        fcCol = findColumn(headers, { ...
            'log2FC', ...
            'log2_FC', ...
            'log2.fc', ...
            'fc'});

        if isempty(idCol)
            idCol = 1;
        end

        if isempty(fcCol)
            error('No log2FC/log2_FC/fc column found in:\n%s', filename);
        end

        ids = raw(2:end, idCol);
        fcValues = raw(2:end, fcCol);

        for j = 1:nMetabolites

            rowIndex = findCellString(ids, selectedIDs{j});

            if isempty(rowIndex)
                warning( ...
                    'Metabolite not found in %s: %s', ...
                    files{f}, ...
                    selectedIDs{j});
                continue;
            end

            result(f,j) = cellToNumber(fcValues{rowIndex});
        end
    end
end

function result = readSelectedReactions( ...
    folder, files, selectedIDs, pathwayName)

    nFiles = numel(files);
    nReactions = numel(selectedIDs);

    result = nan(nFiles, nReactions);

    for f = 1:nFiles

        filename = fullfile(folder, files{f});

        if ~exist(filename, 'file')
            error('File not found:\n%s', filename);
        end

        fprintf('Reading flux-sampling file:\n%s\n', filename);

        [headers, raw] = readExcelRaw(filename);

        idCol = findColumn(headers, { ...
            'Row', ...
            'rxns', ...
            'rxnID', ...
            'ReactionID', ...
            '...1'});

        subsystemCol = findColumn(headers, { ...
            'subSystems', ...
            'subSystem', ...
            'subsys', ...
            'Subsystem'});

        fcCol = findColumn(headers, { ...
            'log2FC', ...
            'log2_FC', ...
            'log2.fc', ...
            'fc'});

        if isempty(idCol)
            idCol = 1;
        end

        if isempty(subsystemCol)
            error('No subsystem column found in:\n%s', filename);
        end

        if isempty(fcCol)
            error('No log2FC/log2_FC/fc column found in:\n%s', filename);
        end

        ids = raw(2:end, idCol);
        subsystems = raw(2:end, subsystemCol);
        fcValues = raw(2:end, fcCol);

        for j = 1:nReactions

            rowIndex = [];

            for r = 1:numel(ids)

                currentID = strtrim(cellToText(ids{r}));
                currentSubsystem = strtrim(cellToText(subsystems{r}));

                if strcmpi(currentID, selectedIDs{j}) && ...
                   strcmpi(currentSubsystem, pathwayName)

                    rowIndex = r;
                    break;
                end
            end

            if isempty(rowIndex)
                warning( ...
                    'Reaction not found in pathway for %s: %s', ...
                    files{f}, ...
                    selectedIDs{j});
                continue;
            end

            result(f,j) = cellToNumber(fcValues{rowIndex});
        end
    end
end

function makePaperDotPlot( ...
    values, xLabels, conditions, plotTitle, outputStub)

    nConditions = size(values,1);
    nFeatures = size(values,2);

    if nConditions ~= 4
        error('Expected four conditions: LIHC, KIRC, Huh7, 769-P.');
    end

    if numel(xLabels) ~= nFeatures
        error('Number of xLabels must match the number of data columns.');
    end

    if numel(conditions) ~= nConditions
        error('Number of conditions must match the number of data rows.');
    end

    yPositions = [4.2, 3.2, 1.2, 0.2];

    colorLimit = 2;
    minMarkerArea = 28;
    maxMarkerArea = 420;
    sizeExponent = 1.10;

    if nFeatures <= 12
        pixelsPerFeature = 72;
        labelFontSize = 12;
        labelRotation = 50;
    elseif nFeatures <= 18
        pixelsPerFeature = 68;
        labelFontSize = 11;
        labelRotation = 55;
    else
        pixelsPerFeature = 72;
        labelFontSize = 10;
        labelRotation = 60;
    end

    figWidth = max(1100, pixelsPerFeature * nFeatures + 420);
    figHeight = 780;

    fig = figure( ...
        'Position', [60, 60, figWidth, figHeight], ...
        'Color', 'w', ...
        'PaperPositionMode', 'auto');

    ax = axes( ...
        'Parent', fig, ...
        'Position', [0.145, 0.31, 0.635, 0.56]);

    hold(ax, 'on');

    colormap(ax, blueGrayRed(256));
    caxis(ax, [-colorLimit colorLimit]);

    for row = 1:nConditions
        for col = 1:nFeatures

            value = values(row,col);

            if isnan(value)
                continue;
            end

            clippedValue = max(min(value, colorLimit), -colorLimit);
            magnitude = abs(clippedValue) / colorLimit;

            markerArea = minMarkerArea + ...
                (maxMarkerArea - minMarkerArea) * ...
                magnitude ^ sizeExponent;

            scatter( ...
                ax, ...
                col, ...
                yPositions(row), ...
                markerArea, ...
                clippedValue, ...
                'filled', ...
                'MarkerEdgeColor', [0.70 0.70 0.70], ...
                'LineWidth', 0.30);
        end
    end

    set(ax, ...
        'XLim', [0.35, nFeatures + 0.65], ...
        'YLim', [-0.35, 4.75], ...
        'XTick', 1:nFeatures, ...
        'XTickLabel', xLabels, ...
        'XTickLabelRotation', labelRotation, ...
        'YTick', [0.2 1.2 3.2 4.2], ...
        'YTickLabel', conditions(end:-1:1), ...
        'TickLabelInterpreter', 'none', ...
        'FontName', 'Arial', ...
        'FontSize', labelFontSize, ...
        'FontWeight', 'bold', ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'LineWidth', 0.8);

    title(ax, plotTitle, ...
        'FontName', 'Arial', ...
        'FontSize', 18, ...
        'FontWeight', 'normal');

    ax.XColor = [0.15 0.15 0.15];
    ax.YColor = [0.15 0.15 0.15];

    bracketX = -1.25;
    textX = -1.72;

    line(ax, [bracketX bracketX], [3.0 4.45], ...
        'Color', [0.15 0.15 0.15], ...
        'LineWidth', 1.3, ...
        'Clipping', 'off');

    line(ax, [bracketX bracketX], [0.0 1.45], ...
        'Color', [0.15 0.15 0.15], ...
        'LineWidth', 1.3, ...
        'Clipping', 'off');

    text(ax, textX, 3.72, 'Patients', ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 15, ...
        'Clipping', 'off');

    text(ax, textX, 0.72, 'Cell lines', ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 15, ...
        'Clipping', 'off');

    cb = colorbar(ax);

    set(cb, ...
        'Position', [0.805, 0.59, 0.024, 0.245], ...
        'Ticks', [-2 -1 0 1 2], ...
        'TickLabels', {'<-2','-1','0','1','>2'}, ...
        'FontName', 'Arial', ...
        'FontSize', 10);

    title(cb, 'log2\_FC', ...
        'Interpreter', 'tex', ...
        'FontName', 'Arial', ...
        'FontSize', 11, ...
        'FontWeight', 'normal');

    legendAx = axes( ...
        'Parent', fig, ...
        'Position', [0.845, 0.535, 0.115, 0.32], ...
        'Visible', 'off');

    hold(legendAx, 'on');

    legendValues = [2 1 0 -1 -2];
    legendLabels = {'>2','1','0','-1','<-2'};
    legendY = 5:-1:1;

    for k = 1:numel(legendValues)

        magnitude = abs(legendValues(k)) / colorLimit;

        markerArea = minMarkerArea + ...
            (maxMarkerArea - minMarkerArea) * ...
            magnitude ^ sizeExponent;

        scatter( ...
            legendAx, ...
            0.32, ...
            legendY(k), ...
            markerArea, ...
            'MarkerEdgeColor', [0.20 0.20 0.20], ...
            'LineWidth', 0.8, ...
            'MarkerFaceColor', 'none');

        text( ...
            legendAx, ...
            0.66, ...
            legendY(k), ...
            legendLabels{k}, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontSize', 10);
    end

    xlim(legendAx, [0 1.25]);
    ylim(legendAx, [0.35 5.65]);

    annotation(fig, 'textbox', ...
        [0.805, 0.425, 0.18, 0.095], ...
        'String', {'Relative to', 'AKR1A1^{high}'}, ...
        'Interpreter', 'tex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 11, ...
        'EdgeColor', 'none');

    annotation(fig, 'textbox', ...
        [0.805, 0.14, 0.18, 0.095], ...
        'String', {'Relative to', 'siSCR'}, ...
        'Interpreter', 'tex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 11, ...
        'EdgeColor', 'none');

    hold(legendAx, 'off');
    hold(ax, 'off');

    svgFile = [outputStub '.svg'];
    pdfFile = [outputStub '.pdf'];

    print(fig, svgFile, '-dsvg');
    print(fig, pdfFile, '-dpdf', '-bestfit');

    fprintf('Created figure:\n%s\n', svgFile);
    fprintf('Created figure:\n%s\n', pdfFile);

    close(fig);
end

function [headers, raw] = readExcelRaw(filename)

    if ~exist(filename, 'file')
        error('Excel file not found:\n%s', filename);
    end

    [~, baseName, extension] = fileparts(filename);

    localCopy = fullfile( ...
        tempdir, ...
        [baseName '_' datestr(now, 'yyyymmdd_HHMMSSFFF') extension]);

    [copyOK, copyMessage] = copyfile(filename, localCopy, 'f');

    if ~copyOK
        error( ...
            ['Could not copy workbook to a local temporary file.\n' ...
             'Source:\n%s\n\nReason:\n%s'], ...
            filename, ...
            copyMessage);
    end

    cleanupObject = onCleanup(@() deleteTemporaryWorkbook(localCopy));

    try
        [~, ~, raw] = xlsread(localCopy, '', '', 'basic');
    catch readError
        error( ...
            ['Could not read workbook:\n%s\n\nMATLAB error:\n%s'], ...
            filename, ...
            readError.message);
    end

    clear cleanupObject;

    if isempty(raw)
        error('Workbook is empty or unreadable:\n%s', filename);
    end

    headers = raw(1,:);
end

function deleteTemporaryWorkbook(filename)

    if exist(filename, 'file')
        try
            delete(filename);
        catch
        end
    end
end

function columnIndex = findColumn(headers, possibleNames)

    columnIndex = [];

    for c = 1:numel(headers)

        currentHeader = strtrim(cellToText(headers{c}));

        for p = 1:numel(possibleNames)

            if strcmpi(currentHeader, possibleNames{p})
                columnIndex = c;
                return;
            end
        end
    end
end

function index = findCellString(cellColumn, target)

    index = [];
    normalizedTarget = normalizeMetaboliteID(target);

    for r = 1:numel(cellColumn)

        currentValue = strtrim(cellToText(cellColumn{r}));

        if strcmpi(currentValue, strtrim(target))
            index = r;
            return;
        end

        normalizedCurrent = normalizeMetaboliteID(currentValue);

        if strcmpi(normalizedCurrent, normalizedTarget)
            index = r;
            return;
        end
    end
end

function normalizedID = normalizeMetaboliteID(value)

    normalizedID = lower(strtrim(cellToText(value)));
    normalizedID = regexprep(normalizedID, '\s+', '');
    normalizedID = regexprep(normalizedID, '\(([a-z])\)$', '[$1]');
    normalizedID = regexprep(normalizedID, '_([a-z])$', '[$1]');
end

function textValue = cellToText(value)

    if isempty(value)
        textValue = '';

    elseif ischar(value)
        textValue = value;

    elseif isstring(value)
        textValue = char(value);

    elseif isnumeric(value) && isscalar(value)
        textValue = num2str(value);

    else
        textValue = '';
    end
end

function numberValue = cellToNumber(value)

    if isnumeric(value) && isscalar(value)
        numberValue = double(value);

    elseif ischar(value)
        numberValue = str2double(strrep(value, ',', '.'));

    elseif isstring(value)
        numberValue = str2double(strrep(char(value), ',', '.'));

    else
        numberValue = NaN;
    end
end

function cmap = blueGrayRed(n)

    if nargin < 1
        n = 256;
    end

    blue = [0.20, 0.31, 0.65];
    gray = [0.78, 0.78, 0.78];
    red  = [0.93, 0.29, 0.20];

    n1 = floor(n / 2);
    n2 = n - n1;

    firstHalf = [ ...
        linspace(blue(1), gray(1), n1)', ...
        linspace(blue(2), gray(2), n1)', ...
        linspace(blue(3), gray(3), n1)' ...
    ];

    secondHalf = [ ...
        linspace(gray(1), red(1), n2)', ...
        linspace(gray(2), red(2), n2)', ...
        linspace(gray(3), red(3), n2)' ...
    ];

    cmap = [firstHalf; secondHalf];
end

function writeSelectedMatrix( ...
    filename, ids, names, conditions, values)

    T = table(ids(:), names(:), ...
        'VariableNames', {'ID', 'Name'});

    for c = 1:numel(conditions)

        validName = makeSafeVariableName(conditions{c});

        if any(strcmp(T.Properties.VariableNames, validName))
            validName = [validName '_data'];
        end

        T.(validName) = values(c,:)';
    end

    if exist(filename, 'file')
        delete(filename);
    end

    writetable(T, filename);
end

function validName = makeSafeVariableName(inputName)

    validName = char(inputName);
    validName = strtrim(validName);
    validName = regexprep(validName, '[^A-Za-z0-9_]', '_');

    if isempty(validName)
        validName = 'Variable';
    end

    if ~isletter(validName(1))
        validName = ['x_' validName];
    end

    if length(validName) > namelengthmax
        validName = validName(1:namelengthmax);
    end
end
