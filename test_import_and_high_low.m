clear;
close all;
clc;

cfg = pipeline_config;

%% Import TCGA annotation file

sampleFile = fullfile( ...
    cfg.rawTCGADir, ...
    'GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt');

if exist(sampleFile, 'file') ~= 2
    error('Missing annotation file: %s', sampleFile);
end

sampleTable = import_TCGA_cancer_type_samples(sampleFile);

disp('TCGA annotation imported.');
disp(size(sampleTable));
disp(sampleTable(1:min(5, height(sampleTable)), :));

%% Make sure output folder exists

if exist(cfg.highLowDir, 'dir') ~= 7
    mkdir(cfg.highLowDir);
end

%% Process TCGA cancer types

cancers = {'KICH', 'KIRC', 'KIRP', 'LIHC'};

highLowResults = struct();

for i = 1:length(cancers)

    cancer = cancers{i};

    fprintf('\nRunning high_low for %s\n', cancer);

    inputFile = fullfile( ...
        cfg.cancerExpressionDir, ...
        [cancer '.txt']);

    if exist(inputFile, 'file') ~= 2
        error('Missing cancer expression file: %s', inputFile);
    end

    %% Read comma-separated expression file

    S = tdfread(inputFile, ',');

    fields = fieldnames(S);

    if length(fields) < 2
        error('No sample columns found in %s', inputFile);
    end

    geneField = fields{1};
    sampleFields = fields(2:end);

    %% Extract genes

    genes = cellstr(S.(geneField));
    genes = strtrim(genes(:));

    nGenes = length(genes);
    nSamples = length(sampleFields);

    %% Build numeric expression matrix

    expressionData = zeros(nGenes, nSamples);

    for k = 1:nSamples

        values = S.(sampleFields{k});

        if ischar(values)
            values = str2double(cellstr(values));
        else
            values = double(values);
        end

        values = values(:);

        if length(values) ~= nGenes
            error( ...
                'Incorrect column length for %s in %s.', ...
                sampleFields{k}, cancer);
        end

        expressionData(:, k) = values;

    end

    %% Run HIGH/LOW selection

    [Table_ON, Table_OFF, colnamesON, colnamesOFF] = high_low( ...
        expressionData, ...
        genes, ...
        sampleFields, ...
        'AKR1A1');

    %% Store results in memory

    highLowResults.(cancer).Table_ON = Table_ON;
    highLowResults.(cancer).Table_OFF = Table_OFF;
    highLowResults.(cancer).Genes = genes;
    highLowResults.(cancer).colnamesON = colnamesON;
    highLowResults.(cancer).colnamesOFF = colnamesOFF;

    %% Output file paths

    fileON = fullfile( ...
        cfg.highLowDir, ...
        ['Table_ON_' cancer '.xlsx']);

    fileOFF = fullfile( ...
        cfg.highLowDir, ...
        ['Table_OFF_' cancer '.xlsx']);

    annoFile = fullfile( ...
        cfg.highLowDir, ...
        [cancer '_anno.txt']);

    %% Delete old files before writing fresh files

    if exist(fileON, 'file') == 2
        delete(fileON);
    end

    if exist(fileOFF, 'file') == 2
        delete(fileOFF);
    end

    if exist(annoFile, 'file') == 2
        delete(annoFile);
    end

    %% Write HIGH and LOW Excel files

    writetable(Table_ON, fileON);
    writetable(Table_OFF, fileOFF);

    %% Write annotation file

    colnamesON = colnamesON(:);
    colnamesOFF = colnamesOFF(:);

    annotationSamples = [colnamesON; colnamesOFF];

    annotationGroups = [ ...
        repmat({'high'}, length(colnamesON), 1); ...
        repmat({'low'}, length(colnamesOFF), 1)];

    annotationTable = table( ...
        annotationSamples, ...
        annotationGroups, ...
        'VariableNames', {'SampleID', 'Group'});

    writetable( ...
        annotationTable, ...
        annoFile, ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'WriteVariableNames', false);

    %% Verify that the Excel files can be reopened

    try
        checkON = readtable(fileON);
    catch ME
        error( ...
            'Table_ON_%s.xlsx was written but cannot be reopened:\n%s', ...
            cancer, ME.message);
    end

    try
        checkOFF = readtable(fileOFF);
    catch ME
        error( ...
            'Table_OFF_%s.xlsx was written but cannot be reopened:\n%s', ...
            cancer, ME.message);
    end

    %% Print results

    fprintf('%s total samples: %d\n', cancer, nSamples);
    fprintf('%s HIGH samples: %d\n', cancer, length(colnamesON));
    fprintf('%s LOW samples: %d\n', cancer, length(colnamesOFF));

    fprintf('Saved ON file: %s\n', fileON);
    fprintf('Saved OFF file: %s\n', fileOFF);
    fprintf('Saved annotation: %s\n', annoFile);

    fprintf('Verified ON file: %d rows, %d columns\n', ...
        height(checkON), width(checkON));

    fprintf('Verified OFF file: %d rows, %d columns\n', ...
        height(checkOFF), width(checkOFF));

end

fprintf('\nHIGH/LOW generation completed successfully.\n');
disp('Results are also available in highLowResults.');