%% driverModel_KICH_ON.m
% Minimal test:
% Build only the KICH ON consensus model.
%
% Original FASTCORMICS settings are preserved.
% Kidney cancer uses RPMI.
% MATLAB R2019b compatible.

clear;
close all;
clc;

feature astheightlimit 2000
solverOK = 1;

cfg = pipeline_config;

%% Temporary test location
% discretizedON_KICH.mat is currently in the project root folder.

cfg.preparedDir = cfg.rootDir;

%% Check required input files

discretizedFile = fullfile( ...
    cfg.preparedDir, ...
    'discretizedON_KICH.mat');

if exist(cfg.modelFile, 'file') ~= 2
    error('Missing model file: %s', cfg.modelFile);
end

if exist(cfg.dicoFile, 'file') ~= 2
    error('Missing dictionary file: %s', cfg.dicoFile);
end

if exist(cfg.mediumFiles.RPMI, 'file') ~= 2
    error('Missing RPMI Excel file: %s', cfg.mediumFiles.RPMI);
end

if exist(discretizedFile, 'file') ~= 2
    error('Missing KICH discretized file: %s', discretizedFile);
end

if exist(cfg.modelsDir, 'dir') ~= 7
    mkdir(cfg.modelsDir);
end

%% Load original metabolic model and dictionary

load(cfg.modelFile, 'consistent_model');
load(cfg.dicoFile, 'dico');

%% Original oxygen setting

consistent_model.lb( ...
    ismember(consistent_model.rxns, 'EX_o2s[e]')) = 0;

%% Original FASTCORMICS settings — unchanged

already_mapped_tag = 0;
consensus_proportion = 0.9;
epsilon = 1e-4;
biomass_rxn = 'biomass_reaction';

%% Original unpenalized systems — unchanged

unpenalizedSystems = { ...
    'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};

unpenalized = consistent_model.rxns( ...
    ismember(consistent_model.subSystems, unpenalizedSystems));

optional_settings.unpenalized = unpenalized;
optional_settings.not_medium_constrained = [];

optional_settings.func = { ...
    'DM_atp_c_', ...
    'biomass_reaction'};

%% Read RPMI medium from Excel

medium_RPMI_7 = readtable( ...
    cfg.mediumFiles.RPMI, ...
    'ReadVariableNames', false);

medium_RPMI_7 = table2array(medium_RPMI_7);

%% Save medium MAT file for later use

mediumMatFile = fullfile( ...
    cfg.modelsDir, ...
    'medium_RPMI_7.mat');

save(mediumMatFile, 'medium_RPMI_7');

%% Assign RPMI medium to KICH

optional_settings.medium = medium_RPMI_7;

%% Load KICH ON discretized data

load( ...
    discretizedFile, ...
    'discretizedon', ...
    'geneSymson', ...
    'colnamesON');

rownames = geneSymson;

%% Validate KICH inputs

if size(discretizedon, 1) ~= numel(rownames)
    error(['KICH ON data mismatch: discretizedon has %d rows, ' ...
           'but geneSymson has %d genes.'], ...
           size(discretizedon, 1), ...
           numel(rownames));
end

if size(discretizedon, 2) ~= numel(colnamesON)
    error(['KICH ON data mismatch: discretizedon has %d columns, ' ...
           'but colnamesON has %d samples.'], ...
           size(discretizedon, 2), ...
           numel(colnamesON));
end

%% Print the exact inputs being used

fprintf('\nBuilding TCGA model: KICH_ON\n');
fprintf('Model: %s\n', cfg.modelFile);
fprintf('Dictionary: %s\n', cfg.dicoFile);
fprintf('Discretized data: %s\n', discretizedFile);
fprintf('RPMI Excel: %s\n', cfg.mediumFiles.RPMI);
fprintf('Genes: %d\n', size(discretizedon, 1));
fprintf('Samples: %d\n', size(discretizedon, 2));
fprintf('Medium entries: %d\n', numel(medium_RPMI_7));
fprintf('Unpenalized reactions: %d\n', numel(unpenalized));

%% Build one KICH ON consensus model

[ContextModel, A_keep] = fastcormics_RNAseq( ...
    consistent_model, ...
    discretizedon, ...
    rownames, ...
    dico, ...
    biomass_rxn, ...
    already_mapped_tag, ...
    consensus_proportion, ...
    epsilon, ...
    optional_settings);

%% Save the model

model_KICH_ON = ContextModel;

models_keep = zeros(numel(consistent_model.rxns), 1);
models_keep(A_keep) = 1;

outputFile = fullfile( ...
    cfg.modelsDir, ...
    'model_KICH_ON.mat');

save( ...
    outputFile, ...
    'model_KICH_ON', ...
    'ContextModel', ...
    'models_keep', ...
    'colnamesON', ...
    'rownames');

fprintf('\nKICH ON model saved:\n%s\n', outputFile);
fprintf('KICH ON test completed.\n');

return