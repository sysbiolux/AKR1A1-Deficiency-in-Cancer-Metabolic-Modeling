%% split_TCGA_high_low.m
%
% Split the original TCGA expression matrices into AKR1A1 HIGH and LOW
% groups for KICH, KIRC, KIRP and LIHC.
%
% Run this script before driverData_shared_clean.m.

clear;
close all;
clc;

%% Relative paths

inputDir = fullfile('data','raw','TCGA_original');
outputDir = fullfile('data','raw','TCGA_ONOFF');

if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

%% Settings

geneName = 'AKR1A1';

cancers = {
    'KICH'
    'KIRC'
    'KIRP'
    'LIHC'
};

% Adapt only these names to the actual original files.
inputFiles = {
    'Table_KICH.xlsx'
    'Table_KIRC.xlsx'
    'Table_KIRP.xlsx'
    'Table_LIHC.xlsx'
};

%% Process each cancer

for cancerIndex = 1:numel(cancers)

    cancer = cancers{cancerIndex};

    fprintf('\n========================================\n');
    fprintf('Splitting %s by %s expression\n',cancer,geneName);
    fprintf('========================================\n');

    inputFile = fullfile( ...
        inputDir, ...
        inputFiles{cancerIndex});

    if exist(inputFile,'file') ~= 2
        error('Original TCGA file not found: %s',inputFile);
    end

    %% Read original unsplit expression table

    T = readtable(inputFile,'PreserveVariableNames',true);

    %% Identify gene column

    variableNames = T.Properties.VariableNames;

    if any(strcmpi(variableNames,'Genes'))
        geneColumn = variableNames{find( ...
            strcmpi(variableNames,'Genes'),1)};
    elseif any(strcmpi(variableNames,'Symbol'))
        geneColumn = variableNames{find( ...
            strcmpi(variableNames,'Symbol'),1)};
    elseif any(strcmpi(variableNames,'Gene'))
        geneColumn = variableNames{find( ...
            strcmpi(variableNames,'Gene'),1)};
    elseif any(strcmpi(variableNames,'Var1'))
        geneColumn = variableNames{find( ...
            strcmpi(variableNames,'Var1'),1)};
    else
        geneColumn = variableNames{1};

        warning( ...
            'Using the first column as the gene column: %s', ...
            geneColumn);
    end

    %% Extract genes

    genes = cellstr(string(T.(geneColumn)));

    %% Extract sample columns

    sampleColumns = variableNames;
    sampleColumns(strcmp(sampleColumns,geneColumn)) = [];

    colnames = sampleColumns(:);

    %% Extract numeric expression data

    data = table2array(T(:,sampleColumns));

    if ~isnumeric(data)
        error( ...
            ['Expression columns in %s could not be converted ' ...
             'into one numeric matrix.'], ...
            inputFile);
    end

    data = double(data);

    %% Perform AKR1A1 HIGH/LOW split

    [Table_ON,Table_OFF,colnamesON,colnamesOFF] = ...
        high_low( ...
            data, ...
            genes, ...
            colnames, ...
            geneName);

    %% Output files

    outputON = fullfile( ...
        outputDir, ...
        ['Table_ON_' cancer '.xlsx']);

    outputOFF = fullfile( ...
        outputDir, ...
        ['Table_OFF_' cancer '.xlsx']);

    writetable(Table_ON,outputON);
    writetable(Table_OFF,outputOFF);

    %% Save the selected names and thresholds-related results

    save( ...
        fullfile(outputDir, ...
        ['high_low_metadata_' cancer '.mat']), ...
        'cancer', ...
        'geneName', ...
        'colnamesON', ...
        'colnamesOFF');

    fprintf('Written: %s\n',outputON);
    fprintf('Written: %s\n',outputOFF);

end

disp('TCGA HIGH/LOW splitting finished');