clear; close all; clc;

cfg = pipeline_config;

sampleTable = import_TCGA_cancer_type_samples();

cancers = {'KICH','KIRC','KIRP','LIHC'};

for i = 1:length(cancers)

    cancer = cancers{i};
    disp(['Processing ' cancer]);

    inputFile = fullfile(cfg.cancerExpressionDir, [cancer '.txt']);

    if exist(inputFile,'file') ~= 2
        warning(['Missing file: ' inputFile]);
        continue
    end

    T = readtable(inputFile, ...
        'FileType','text', ...
        'Delimiter','\t', ...
        'ReadVariableNames',true);

    [Table_ON, Table_OFF, colnamesON, colnamesOFF] = high_low(T,'AKR1A1');

    writetable(Table_ON, ...
        fullfile(cfg.highLowDir, ['Table_ON_' cancer '.xlsx']));

    writetable(Table_OFF, ...
        fullfile(cfg.highLowDir, ['Table_OFF_' cancer '.xlsx']));

    annoSamples = [colnamesON, colnamesOFF];
    annoGroups = [repmat({'high'},1,length(colnamesON)), ...
                  repmat({'low'},1,length(colnamesOFF))];

    annoTable = cell2table([annoSamples; annoGroups]);

    writetable(annoTable, ...
        fullfile(cfg.highLowDir, [cancer '_anno.txt']), ...
        'FileType','text', ...
        'Delimiter','\t', ...
        'WriteVariableNames',false);

    disp(['Finished ' cancer]);

end