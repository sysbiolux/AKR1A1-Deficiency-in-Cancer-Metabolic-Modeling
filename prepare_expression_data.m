clear; close all; clc;

cfg = pipeline_config;

if strcmpi(cfg.runMode,'TCGA') || strcmpi(cfg.runMode,'ALL')
    prepare_TCGA(cfg);
end

if strcmpi(cfg.runMode,'CELL_LINES') || strcmpi(cfg.runMode,'ALL')
    prepare_cell_lines(cfg);
end

disp('Step 02 finished');

function prepare_TCGA(cfg)

for c = 1:numel(cfg.cancerTypes)

    cancer = cfg.cancerTypes{c};

    onFile = fullfile(cfg.highLowDir, ['Table_ON_' cancer '.xlsx']);
    offFile = fullfile(cfg.highLowDir, ['Table_OFF_' cancer '.xlsx']);

    if exist(onFile,'file') ~= 2 || exist(offFile,'file') ~= 2
        warning(['Missing ON/OFF files for ' cancer]);
        continue
    end

    T_on = readtable(onFile);
    T_off = readtable(offFile);

    [discretizedon, geneSymson, colnamesON] = prepare_table(T_on);
    [discretizedoff, geneSymsoff, colnamesOFF] = prepare_table(T_off);

    save(fullfile(cfg.preparedDir, ['discretizedON_' cancer '.mat']), ...
        'discretizedon','geneSymson','colnamesON');

    save(fullfile(cfg.preparedDir, ['discretizedOFF_' cancer '.mat']), ...
        'discretizedoff','geneSymsoff','colnamesOFF');

    disp(['Prepared ' cancer]);
end

end

function prepare_cell_lines(cfg)

keys = fieldnames(cfg.cellLines);

for k = 1:numel(keys)

    key = keys{k};
    info = cfg.cellLines.(key);

    if exist(info.expression,'file') ~= 2
        warning(['Missing cell-line file: ' info.expression]);
        continue
    end

    T = readtable(info.expression, ...
        'FileType','text', ...
        'Delimiter','\t', ...
        'ReadVariableNames',true);

    [discretized, geneSyms, colnames] = prepare_table(T);

    outFile = fullfile(cfg.preparedDir, ['discretized_' info.name '.mat']);
    save(outFile, 'discretized','geneSyms','colnames','info');

    disp(['Prepared ' info.name]);
end

end

function [discretized, geneSyms, colnames] = prepare_table(T)

varNames = T.Properties.VariableNames;

geneCol = find(strcmpi(varNames,'Var1') | strcmpi(varNames,'Genes') | strcmpi(varNames,'Symbol'), 1);

if isempty(geneCol)
    geneCol = 1;
end

geneSyms = cellstr(string(T{:,geneCol}));

sampleCols = setdiff(1:width(T), geneCol);
colnames = T.Properties.VariableNames(sampleCols);

data = T{:,sampleCols};

if ~isnumeric(data)
    data = str2double(string(data));
end

data(isnan(data)) = 0;

keep = ~strcmpi(geneSyms,'NA') & ~strcmpi(geneSyms,'');
geneSyms = geneSyms(keep);
data = data(keep,:);

if exist('discretize_FPKM','file') == 2
    discretized = discretize_FPKM(data, colnames, 1);
else
    error('discretize_FPKM was not found on the MATLAB path');
end

end