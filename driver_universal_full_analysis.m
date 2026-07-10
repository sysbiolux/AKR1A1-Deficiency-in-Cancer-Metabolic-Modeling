%% DRIVER_UNIVERSAL_FULL_ANALYSIS
%
% Run the complete universal analysis for all models in one dataset.
%
% MATLAB version:
%   R2019b

clear;
close all;
clc;

feature astheightlimit 2000

%% ============================================================
% USER SETTING
% =============================================================

% Available options:
%
%   'TCGA_ONOFF'
%   'CELLLINES'

datasetType = 'TCGA_ONOFF';

%% ============================================================
% LOAD DATASET CONFIGURATION
% =============================================================

switch upper(datasetType)

    case 'TCGA_ONOFF'

        config = config_TCGA_ONOFF();

    case 'CELLLINES'

        config = config_CellLines();

    otherwise

        error( ...
            ['Unknown datasetType: %s\n' ...
             'Use TCGA_ONOFF or CELLLINES.'], ...
            datasetType);

end

%% ============================================================
% SELECT COMPLETE ANALYSIS
% =============================================================

config.runSingleGeneDeletion = true;
config.runDrugDeletion       = true;
config.runFVA                = true;
config.runFVASimilarity      = true;

config.saveIndividualGeneDeletionFiles = true;
config.saveIndividualDrugFiles         = true;
config.saveCompleteResults             = true;

%% ============================================================
% OUTPUT LOCATION
% =============================================================

config.outputFolder = fullfile( ...
    pwd, ...
    'results', ...
    config.datasetName);

%% ============================================================
% DISPLAY SETTINGS
% =============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('Complete universal analysis\n');
fprintf('============================================\n');

fprintf('Dataset: %s\n', config.datasetName);
fprintf('Number of models: %d\n', numel(config.models));

fprintf('Models:\n');

for modelIndex = 1:numel(config.modelNames)

    fprintf( ...
        '  %d. %s\n', ...
        modelIndex, ...
        config.modelNames{modelIndex});

end

fprintf('\nOutput folder:\n');
fprintf('  %s\n', config.outputFolder);

%% ============================================================
% RUN ANALYSIS
% =============================================================

results = universal_model_analysis(config);

fprintf('\n');
fprintf('============================================\n');
fprintf('COMPLETE ANALYSIS FINISHED\n');
fprintf('============================================\n');

fprintf('Results saved in:\n');
fprintf('  %s\n', config.outputFolder);