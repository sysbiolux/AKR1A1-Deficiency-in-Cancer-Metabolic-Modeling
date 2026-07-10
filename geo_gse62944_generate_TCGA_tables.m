clear; close all; clc;

cfg = pipeline_config;

sampleFile = fullfile(cfg.rawTCGADir, cfg.tcga.sampleAnnotationFile);
tumorFile = fullfile(cfg.rawTCGADir, cfg.tcga.tumorExpressionFile);

sampleTable = import_TCGA_cancer_type_samples(sampleFile);
sampleTable.SampleIDClean = clean_ids(sampleTable.SampleID);

fid = fopen(tumorFile,'r');
headerLine = fgetl(fid);
fclose(fid);

headers = strsplit(headerLine, sprintf('\t'));
geneHeader = headers{1};
exprSampleNames = headers(2:end);
exprSampleClean = clean_ids(exprSampleNames);

for c = 1:numel(cfg.cancerTypes)

    cancer = cfg.cancerTypes{c};
    disp(['Streaming cancer file: ' cancer]);

    cancerMask = strcmp(cellstr(sampleTable.CancerType), cancer);
    cancerSamples = sampleTable.SampleIDClean(cancerMask);

    matchMask = ismember(exprSampleClean, cancerSamples);
    selectedCols = find(matchMask) + 1;

    if isempty(selectedCols)
        warning(['No samples found for ' cancer]);
        continue
    end

    outFile = fullfile(cfg.cancerExpressionDir, [cancer '.txt']);

    fin = fopen(tumorFile,'r');
    fout = fopen(outFile,'w');

    header = strsplit(fgetl(fin), sprintf('\t'));
    fprintf(fout, 'Genes');
    for j = 1:numel(selectedCols)
        fprintf(fout, '\t%s', header{selectedCols(j)});
    end
    fprintf(fout, '\n');

    while ~feof(fin)
        line = fgetl(fin);
        if ~ischar(line)
            continue
        end

        parts = strsplit(line, sprintf('\t'));

        fprintf(fout, '%s', parts{1});

        for j = 1:numel(selectedCols)
            col = selectedCols(j);
            if col <= numel(parts)
                fprintf(fout, '\t%s', parts{col});
            else
                fprintf(fout, '\t0');
            end
        end

        fprintf(fout, '\n');
    end

    fclose(fin);
    fclose(fout);

    disp(['Created ' outFile]);

    T = readtable(outFile, 'FileType','text', 'Delimiter','\t', 'ReadVariableNames',true);

    [Table_ON, Table_OFF, annoTable] = high_low_table(T, cfg.highLow.gene, ...
        cfg.highLow.highPercentile, cfg.highLow.lowPercentile);

    writetable(Table_ON, fullfile(cfg.highLowDir, ['Table_ON_' cancer '.xlsx']));
    writetable(Table_OFF, fullfile(cfg.highLowDir, ['Table_OFF_' cancer '.xlsx']));
    writetable(annoTable, fullfile(cfg.highLowDir, [cancer '_anno.txt']), ...
        'FileType','text','Delimiter','\t','WriteVariableNames',false);

    disp(['Created high-low files for ' cancer]);
end

disp('Step 01 finished');

function sampleTable = import_TCGA_cancer_type_samples(filename)

opts = delimitedTextImportOptions('NumVariables', 2);
opts.DataLines = [1, Inf];
opts.Delimiter = '\t';
opts.VariableNames = {'SampleID','CancerType'};
opts.VariableTypes = {'char','categorical'};
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';
opts = setvaropts(opts, 'SampleID', 'WhitespaceRule', 'preserve');
opts = setvaropts(opts, {'SampleID','CancerType'}, 'EmptyFieldRule', 'auto');

sampleTable = readtable(filename, opts);

end

function cleanIDs = clean_ids(ids)

ids = cellstr(ids);
cleanIDs = cell(size(ids));

for i = 1:numel(ids)
    x = ids{i};
    x = strrep(x, '-', '');
    x = strrep(x, '.', '');
    x = strrep(x, '_', '');
    cleanIDs{i} = upper(x);
end

end

function [Table_ON, Table_OFF, annoTable] = high_low_table(T, geneName, highPct, lowPct)

genes = cellstr(T.Genes);
geneIdx = find(strcmpi(genes, geneName));

if isempty(geneIdx)
    geneIdx = find(contains(upper(genes), upper(geneName)));
end

if isempty(geneIdx)
    error(['Gene not found: ' geneName]);
end

geneIdx = geneIdx(1);

sampleNames = T.Properties.VariableNames(2:end);
data = table2array(T(:,2:end));

akr = data(geneIdx,:);

selectON = akr >= prctile(akr, highPct);
selectOFF = akr <= prctile(akr, lowPct);

Table_ON = T(:, [{'Genes'}, sampleNames(selectON)]);
Table_OFF = T(:, [{'Genes'}, sampleNames(selectOFF)]);

Table_ON.Properties.VariableNames{1} = 'Var1';
Table_OFF.Properties.VariableNames{1} = 'Var1';

annoSamples = [sampleNames(selectON), sampleNames(selectOFF)];
annoGroups = [repmat({'high'},1,sum(selectON)), repmat({'low'},1,sum(selectOFF))];

annoTable = cell2table([annoSamples; annoGroups]);

end