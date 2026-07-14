%% driverData_shared_clean.m
%
% Shared data driver for:
%
% 1) Standard cell-line data:
%       Huh7 and 769-P
%
% 2) TCGA AKR1A1 ON/OFF data:
%       KICH, KIRC, KIRP and LIHC
%
% Before running the TCGA branch, run:
%
%   split_TCGA_high_low
%
% MATLAB R2019b compatible.

clear;
close all;
clc;

solverOK = 1;

%% ============================================================
% EXECUTION MODE
%% ============================================================

% Available options:
%   'CELL_LINES'
%   'TCGA'
%   'ALL'

runMode = 'ALL';

%% ============================================================
% RELATIVE PATHS
%% ============================================================

dataRawDir = fullfile('data','raw');
preparedDir = fullfile('data','prepared');
imagesDir = 'images';

huh7Dir = fullfile(dataRawDir,'HuH7');
p769Dir = fullfile(dataRawDir,'769-p');
highLowDir = fullfile(dataRawDir,'TCGA_ONOFF');

geneLengthFile = fullfile( ...
    dataRawDir, ...
    'geneLengths_170423.xlsx');

%% Create output directories

if exist(preparedDir,'dir') ~= 7
    mkdir(preparedDir);
end

if exist(imagesDir,'dir') ~= 7
    mkdir(imagesDir);
end

%% Validate execution mode

validModes = {'CELL_LINES','TCGA','ALL'};

if ~any(strcmpi(runMode,validModes))
    error('Invalid runMode: %s',runMode);
end

%% Check required discretization function

if exist('discretize_FPKM','file') ~= 2
    error('discretize_FPKM.m was not found on the MATLAB path.');
end

%% ============================================================
% TCGA AKR1A1 ON/OFF DATA
%% ============================================================

if strcmpi(runMode,'TCGA') || strcmpi(runMode,'ALL')

    if exist('read_excel_data','file') ~= 2
        error('read_excel_data.m was not found on the MATLAB path.');
    end

    filesOFF = {
        'Table_OFF_KICH.xlsx'
        'Table_OFF_KIRC.xlsx'
        'Table_OFF_KIRP.xlsx'
        'Table_OFF_LIHC.xlsx'
    };

    filesON = {
        'Table_ON_KICH.xlsx'
        'Table_ON_KIRC.xlsx'
        'Table_ON_KIRP.xlsx'
        'Table_ON_LIHC.xlsx'
    };

    cancers = {
        'KICH'
        'KIRC'
        'KIRP'
        'LIHC'
    };

    for cancerIndex = 1:numel(cancers)

        cancer = cancers{cancerIndex};

        fprintf('\n========================================\n');
        fprintf('Preparing TCGA %s\n',cancer);
        fprintf('========================================\n');

        file_off = fullfile( ...
            highLowDir, ...
            filesOFF{cancerIndex});

        file_on = fullfile( ...
            highLowDir, ...
            filesON{cancerIndex});

        %% Check split files

        if exist(file_off,'file') ~= 2
            error(['TCGA OFF file not found: %s\n' ...
                   'Run split_TCGA_high_low first.'], ...
                   file_off);
        end

        if exist(file_on,'file') ~= 2
            error(['TCGA ON file not found: %s\n' ...
                   'Run split_TCGA_high_low first.'], ...
                   file_on);
        end

        %% Read split files as structures

        temp_off = read_excel_data(file_off);
        temp_on = read_excel_data(file_on);

        %% Identify gene-symbol fields

        symbolFieldOff = get_symbol_field(temp_off);
        symbolFieldOn = get_symbol_field(temp_on);

        geneSymsoff = structure_field_to_cellstr( ...
            temp_off.(symbolFieldOff));

        geneSymson = structure_field_to_cellstr( ...
            temp_on.(symbolFieldOn));

        %% Identify expression fields

        fieldsOFF = fieldnames(temp_off);
        fieldsON = fieldnames(temp_on);

        colnamesOFF = fieldsOFF(~strcmp(fieldsOFF,symbolFieldOff));
        colnamesON = fieldsON(~strcmp(fieldsON,symbolFieldOn));

        %% Build OFF expression matrix

        data_valuesoff = [];

        for sampleIndex = 1:numel(colnamesOFF)

            currentValues = temp_off.(colnamesOFF{sampleIndex});
            currentValues = double(currentValues);

            if size(currentValues,2) > 1
                currentValues = currentValues(:,1);
            end

            data_valuesoff = [data_valuesoff,currentValues]; %#ok<AGROW>
        end

        %% Build ON expression matrix

        data_valueson = [];

        for sampleIndex = 1:numel(colnamesON)

            currentValues = temp_on.(colnamesON{sampleIndex});
            currentValues = double(currentValues);

            if size(currentValues,2) > 1
                currentValues = currentValues(:,1);
            end

            data_valueson = [data_valueson,currentValues]; %#ok<AGROW>
        end

        %% Check dimensions

        if size(data_valuesoff,1) ~= numel(geneSymsoff)
            error('%s OFF gene and expression row counts differ.',cancer);
        end

        if size(data_valueson,1) ~= numel(geneSymson)
            error('%s ON gene and expression row counts differ.',cancer);
        end

        if numel(geneSymsoff) ~= numel(geneSymson)

            warning('%s ON and OFF contain different numbers of genes.', ...
                cancer);

        elseif ~isequal(geneSymsoff,geneSymson)

            warning('%s ON and OFF gene ordering is not identical.', ...
                cancer);
        end

        %% Identify entries without annotated symbols

        isPseudoOff = find(strcmpi(geneSymsoff,'NA'));
        isPseudoOn = find(strcmpi(geneSymson,'NA'));

        fprintf('Pseudo genes OFF: %d\n',numel(isPseudoOff));
        fprintf('Pseudo genes ON: %d\n',numel(isPseudoOn));

        %% Discretize ON/OFF expression

        discretizedon = discretize_FPKM( ...
            data_valueson, ...
            colnamesON, ...
            1);

        discretizedoff = discretize_FPKM( ...
            data_valuesoff, ...
            colnamesOFF, ...
            1);

        %% Save separate files used by model reconstruction

        save( ...
            fullfile(preparedDir, ...
            ['discretizedON_' cancer '.mat']), ...
            'discretizedon', ...
            'geneSymson', ...
            'colnamesON');

        save( ...
            fullfile(preparedDir, ...
            ['discretizedOFF_' cancer '.mat']), ...
            'discretizedoff', ...
            'geneSymsoff', ...
            'colnamesOFF');

        %% ON discretization statistics

        on = sum(discretizedon == 1);
        nd = sum(discretizedon == 0);
        off = sum(discretizedon == -1);

        fprintf('ON statistics for %s\n',cancer);
        disp(on);
        disp(nd);
        disp(off);
        disp(on + nd + off);

        %% OFF discretization statistics

        on = sum(discretizedoff == 1);
        nd = sum(discretizedoff == 0);
        off = sum(discretizedoff == -1);

        fprintf('OFF statistics for %s\n',cancer);
        disp(on);
        disp(nd);
        disp(off);
        disp(on + nd + off);

        %% Save complete cancer-specific data workspace

        save( ...
            fullfile(preparedDir, ...
            ['driverDataOFFON_' cancer '.mat']), ...
            'data_valueson', ...
            'data_valuesoff', ...
            'geneSymson', ...
            'geneSymsoff', ...
            'colnamesON', ...
            'colnamesOFF', ...
            'discretizedon', ...
            'discretizedoff');

    end
end

%% ============================================================
% STANDARD CELL-LINE DATA
%% ============================================================

if strcmpi(runMode,'CELL_LINES') || strcmpi(runMode,'ALL')

    fprintf('\n========================================\n');
    fprintf('Preparing Huh7 and 769-P data\n');
    fprintf('========================================\n');

    %% Define input files

    huh7File = fullfile( ...
        huh7Dir, ...
        'Huh7_project_1_normalizedCountsWithAnnotations.txt');

    p769File = fullfile( ...
        p769Dir, ...
        '769P_project_1_normalizedCountsWithAnnotations.txt');

    %% Check input files

    if exist(huh7File,'file') ~= 2
        error('Huh7 input file not found: %s',huh7File);
    end

    if exist(p769File,'file') ~= 2
        error('769-P input file not found: %s',p769File);
    end

    if exist(geneLengthFile,'file') ~= 2

        alternativeGeneLengthFile = 'geneLengths_170423.xlsx';

        if exist(alternativeGeneLengthFile,'file') == 2
            geneLengthFile = alternativeGeneLengthFile;
        else
            error('Gene-length file not found: %s',geneLengthFile);
        end
    end

    %% Read normalized expression files

    temp = tdfread(huh7File);
    temp2 = tdfread(p769File);

    %% Extract sample names using the original field positions

    tempFields = fieldnames(temp);
    colnamesH = tempFields(8:18);

    tempFields = fieldnames(temp2);
    colnames7 = tempFields(8:18);

    %% Separate pseudogene indices

    isPseudoH = find(strcmpi(cellstr(temp.Symbol),'NA'));
    isPseudo7 = find(strcmpi(cellstr(temp2.Symbol),'NA'));

    fprintf('Pseudo genes Huh7: %d\n',numel(isPseudoH));
    fprintf('Pseudo genes 769-P: %d\n',numel(isPseudo7));

    %% Original Huh7 sample ordering

    data_h = [ ...
        temp.Huh7_siSCR_1, ...
        temp.Huh7_siSCR_2, ...
        temp.Huh7_siSCR_3, ...
        temp.Huh7_siSCR_4];

    data_h = [ ...
        data_h, ...
        temp.Huh7_siAKR1A1_1_2, ...
        temp.Huh7_siAKR1A1_1_3, ...
        temp.Huh7_siAKR1A1_1_4];

    data_h = [ ...
        data_h, ...
        temp.Huh7_siAKR1A1_2_1, ...
        temp.Huh7_siAKR1A1_2_2, ...
        temp.Huh7_siAKR1A1_2_3, ...
        temp.Huh7_siAKR1A1_2_4];

    %% Original 769-P sample ordering

    data_r = [ ...
        temp2.x7690x2DP_siSCR_1, ...
        temp2.x7690x2DP_siSCR_2, ...
        temp2.x7690x2DP_siSCR_3, ...
        temp2.x7690x2DP_siSCR_4];

    data_r = [ ...
        data_r, ...
        temp2.x7690x2DP_siAKR1A1_1_1, ...
        temp2.x7690x2DP_siAKR1A1_1_2, ...
        temp2.x7690x2DP_siAKR1A1_1_3];

    data_r = [ ...
        data_r, ...
        temp2.x7690x2DP_siAKR1A1_2_1, ...
        temp2.x7690x2DP_siAKR1A1_2_2, ...
        temp2.x7690x2DP_siAKR1A1_2_3, ...
        temp2.x7690x2DP_siAKR1A1_2_4];

    %% Gene identifiers and symbols

    genes_id_H = cellstr(temp.Alias);
    genes_id_7 = cellstr(temp2.Alias);

    geneSyms_orig_H = cellstr(temp.Symbol);
    geneSyms_orig_7 = cellstr(temp2.Symbol);

    %% Read gene lengths

    [NUM,~,~] = xlsread(geneLengthFile);
    lengths_orig = NUM(:,1);

    if numel(lengths_orig) ~= size(data_h,1)
        error(['Gene-length rows do not match the Huh7 ' ...
               'expression rows.']);
    end

    if numel(lengths_orig) ~= size(data_r,1)
        error(['Gene-length rows do not match the 769-P ' ...
               'expression rows.']);
    end

    %% Process Huh7 and 769-P

    for datasetIndex = 1:2

        if datasetIndex == 1

            data = data_h;
            geneSyms = geneSyms_orig_H;
            genes_id = genes_id_H;
            colnamesCurrent = colnamesH;
            isPseudo = isPseudoH;

            outputName = 'discretizedH.mat';
            datasetLabel = 'Huh7';

        else

            data = data_r;
            geneSyms = geneSyms_orig_7;
            genes_id = genes_id_7;
            colnamesCurrent = colnames7;
            isPseudo = isPseudo7;

            outputName = 'discretized7.mat';
            datasetLabel = '769-P';

        end

        fprintf('\nPreparing %s\n',datasetLabel);

        %% Return transformed expression to original scale

        data = 2.^data - 1;

        %% Set entries without annotated gene symbols to zero

        data(isPseudo,:) = 0;

        lengths = lengths_orig;

        %% Remove genes without gene-length information

        remove = isnan(lengths);

        if any(remove)

            data(remove,:) = [];
            lengths(remove) = [];
            genes_id(remove) = [];
            geneSyms(remove) = [];
        end

        if any(lengths == 0)
            error('%s contains genes with zero gene length.',datasetLabel);
        end

        %% TPM calculation

        tempTPM = data ./ repmat( ...
            lengths, ...
            1, ...
            size(data,2));

        TPM = 1e6 .* tempTPM ./ repmat( ...
            nansum(tempTPM), ...
            size(data,1), ...
            1);

        %% FPKM calculation

        tempFPKM = 1e6 .* data ./ repmat( ...
            nansum(data), ...
            size(data,1), ...
            1);

        FPKM = tempFPKM ./ repmat( ...
            lengths, ...
            1, ...
            size(data,2));

        %% Original inspection object

        data3 = TPM;
        data3(data3 == 0) = NaN;
        data3 = log2(data3);

        %% Discretize TPM data

        discretized = discretize_FPKM( ...
            TPM, ...
            colnamesCurrent, ...
            1);

        %% Save original variable names

        if datasetIndex == 1

            colnamesH = colnamesCurrent;

            save( ...
                fullfile(preparedDir,outputName), ...
                'TPM', ...
                'FPKM', ...
                'geneSyms', ...
                'genes_id', ...
                'colnamesH', ...
                'discretized');

        else

            colnames7 = colnamesCurrent;

            save( ...
                fullfile(preparedDir,outputName), ...
                'TPM', ...
                'FPKM', ...
                'geneSyms', ...
                'genes_id', ...
                'colnames7', ...
                'discretized');

        end

        %% Discretization statistics

        on = sum(discretized == 1);
        nd = sum(discretized == 0);
        off = sum(discretized == -1);

        fprintf('Discretization statistics for %s\n',datasetLabel);
        disp(on);
        disp(nd);
        disp(off);
        disp(on + nd + off);

    end

    %% Save standard cell-line workspace

    save( ...
        fullfile(preparedDir,'driverData_cell_lines.mat'), ...
        'data_h', ...
        'data_r', ...
        'colnamesH', ...
        'colnames7', ...
        'geneSyms_orig_H', ...
        'geneSyms_orig_7');

end

%% Save final workspace

save(fullfile( ...
    preparedDir, ...
    'driverData_shared_clean_workspace.mat'));

disp('driverData_shared_clean finished');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function symbolField = get_symbol_field(S)
%GET_SYMBOL_FIELD Identify the gene-symbol field in a structure.

fields = fieldnames(S);

preferredFields = {
    'Genes'
    'Gene'
    'Symbol'
    'GeneSymbol'
    'Var1'
};

symbolField = '';

for i = 1:numel(preferredFields)

    match = strcmpi(fields,preferredFields{i});

    if any(match)
        symbolField = fields{find(match,1,'first')};
        return
    end
end

symbolField = fields{1};

warning(['No recognized gene-symbol field was found. ' ...
         'Using the first field: %s'], ...
         symbolField);

end


function values = structure_field_to_cellstr(rawValues)
%STRUCTURE_FIELD_TO_CELLSTR Convert a tdfread field to trimmed cell strings.

if ischar(rawValues)
    values = cellstr(rawValues);
elseif iscell(rawValues)
    values = rawValues;
elseif isstring(rawValues)
    values = cellstr(rawValues);
elseif iscategorical(rawValues)
    values = cellstr(rawValues);
else
    values = cellstr(string(rawValues));
end

values = values(:);

for i = 1:numel(values)
    values{i} = strtrim(char(values{i}));
end

end