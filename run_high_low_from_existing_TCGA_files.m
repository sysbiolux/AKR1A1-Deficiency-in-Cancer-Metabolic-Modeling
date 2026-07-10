clear; close all; clc;

cfg = pipeline_config;

sampleFile = fullfile(cfg.rawTCGADir, ...
    'GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt');

sampleTable = import_TCGA_cancer_type_samples(sampleFile);

writetable(sampleTable, ...
    fullfile(cfg.resultsDir,'TCGA_24_CancerType_Samples_imported.csv'));

cancers = {'KICH','KIRC','KIRP','LIHC'};

for i = 1:length(cancers)

    cancer = cancers{i};
    disp(['Running high_low for ' cancer]);

    inputFile = fullfile(cfg.cancerExpressionDir,[cancer '.txt']);

    S = tdfread(inputFile);

    fields = fieldnames(S);
    geneField = fields{1};
    sampleFields = fields(2:end);

    genes = cellstr(S.(geneField));

    data = zeros(length(genes), length(sampleFields));
    for k = 1:length(sampleFields)
        data(:,k) = double(S.(sampleFields{k}));
    end

    [Table_ON, Table_OFF, colnamesON, colnamesOFF] = ...
        high_low(data, genes, sampleFields, 'AKR1A1');

    writetable(Table_ON, fullfile(cfg.highLowDir, ['Table_ON_' cancer '.xlsx']));
    writetable(Table_OFF, fullfile(cfg.highLowDir, ['Table_OFF_' cancer '.xlsx']));

    annoSamples = [colnamesON, colnamesOFF];
    annoGroups = [repmat({'high'},1,length(colnamesON)), ...
                  repmat({'low'},1,length(colnamesOFF))];

    annoTable = cell2table([annoSamples; annoGroups]);

    writetable(annoTable, fullfile(cfg.highLowDir, [cancer '_anno.txt']), ...
        'Delimiter','\t', 'WriteVariableNames',false);
    

    disp(['Finished ' cancer]);
end