%% DRIVER_UNIVERSAL_TEST
%
% Test the universal metabolic-model analysis pipeline in stages.
%
% MATLAB version:
%   R2019b
%
% This driver allows testing:
%
%   Stage 1: FVA with one model
%   Stage 2: Drug deletion with one model
%   Stage 3: Single-gene deletion with one model
%   Stage 4: All single-model analyses
%   Stage 5: FVA similarity with two models
%
% The scientific algorithms are implemented in:
%
%   universal_model_analysis.m
%
% This driver only selects the dataset, models, and analyses.

clear;
close all;
clc;

feature astheightlimit 2000

fprintf('\n');
fprintf('============================================\n');
fprintf('Universal metabolic analysis test driver\n');
fprintf('MATLAB version: %s\n', version);
fprintf('============================================\n');

%% ============================================================
% USER SETTINGS
% =============================================================

%% Select the dataset
%
% Available options:
%
%   'TCGA_ONOFF'
%   'CELLLINES'

datasetType = 'TCGA_ONOFF';

%% Select the test stage
%
%   1 = One model: FVA only
%   2 = One model: drug deletion only
%   3 = One model: single-gene deletion only
%   4 = One model: FVA + drug deletion + gene deletion
%   5 = Two models: FVA + FVA similarity

testStage = 1;

%% Select the model or models
%
% For stages 1 to 4, only the first value is used.
%
% TCGA configuration:
%
%   1 = KIRC_ON
%   2 = KIRC_OFF
%   3 = KIRP_ON
%   4 = KIRP_OFF
%   5 = KICH_ON
%   6 = KICH_OFF
%   7 = LIHC_ON
%   8 = LIHC_OFF
%
% Cell-line configuration:
%
%   1 = 769P_ctrl
%   2 = 769P_sc1
%   3 = 769P_sc2
%   4 = 769P_sc12
%   5 = Huh7_ctrl
%   6 = Huh7_sc1
%   7 = Huh7_sc2
%   8 = Huh7_sc12
%
% Examples:
%
%   One-model test:
%       testModelIndices = 1;
%
%   Two-model similarity test:
%       testModelIndices = [1, 2];

testModelIndices = [1, 2];

%% ============================================================
% LOAD DATASET CONFIGURATION
% =============================================================

switch upper(datasetType)

    case 'TCGA_ONOFF'

        fullConfig = config_TCGA_ONOFF();

    case 'CELLLINES'

        fullConfig = config_CellLines();

    otherwise

        error( ...
            ['Unknown datasetType: %s\n' ...
             'Use TCGA_ONOFF or CELLLINES.'], ...
            datasetType);

end

%% ============================================================
% CHECK TEST-STAGE VALUE
% =============================================================

if ~ismember(testStage, 1:5)

    error( ...
        'testStage must be an integer between 1 and 5.');

end

%% ============================================================
% DETERMINE NUMBER OF REQUIRED MODELS
% =============================================================

if testStage <= 4

    numberOfTestModels = 1;

else

    numberOfTestModels = 2;

end

%% ============================================================
% CHECK MODEL SELECTION
% =============================================================

if numel(testModelIndices) < numberOfTestModels

    error( ...
        ['Test stage %d requires %d model index or indices, ' ...
         'but only %d were provided.'], ...
        testStage, ...
        numberOfTestModels, ...
        numel(testModelIndices));

end

testModelIndices = ...
    testModelIndices(1:numberOfTestModels);

numberOfAvailableModels = ...
    numel(fullConfig.models);

if any(testModelIndices < 1) || ...
        any(testModelIndices > numberOfAvailableModels)

    error( ...
        ['The selected model indices must be between 1 and %d.'], ...
        numberOfAvailableModels);

end

if numel(unique(testModelIndices)) ~= ...
        numel(testModelIndices)

    error('The test model indices must be different.');

end

%% ============================================================
% CREATE THE TEST CONFIGURATION
% =============================================================

config = fullConfig;

config.models = ...
    fullConfig.models(testModelIndices);

config.modelNames = ...
    fullConfig.modelNames(testModelIndices);

if isfield(fullConfig, 'groupNames')

    config.groupNames = ...
        fullConfig.groupNames(testModelIndices);

end

if isfield(fullConfig, 'conditionNames')

    config.conditionNames = ...
        fullConfig.conditionNames(testModelIndices);

end

%% ============================================================
% SELECT ANALYSES ACCORDING TO TEST STAGE
% =============================================================

switch testStage

    case 1

        stageName = 'ONE_MODEL_FVA';

        config.runSingleGeneDeletion = false;
        config.runDrugDeletion       = false;
        config.runFVA                = true;
        config.runFVASimilarity      = false;

    case 2

        stageName = 'ONE_MODEL_DRUG_DELETION';

        config.runSingleGeneDeletion = false;
        config.runDrugDeletion       = true;
        config.runFVA                = false;
        config.runFVASimilarity      = false;

    case 3

        stageName = 'ONE_MODEL_GENE_DELETION';

        config.runSingleGeneDeletion = true;
        config.runDrugDeletion       = false;
        config.runFVA                = false;
        config.runFVASimilarity      = false;

    case 4

        stageName = 'ONE_MODEL_COMPLETE';

        config.runSingleGeneDeletion = true;
        config.runDrugDeletion       = true;
        config.runFVA                = true;
        config.runFVASimilarity      = false;

    case 5

        stageName = 'TWO_MODEL_FVA_SIMILARITY';

        config.runSingleGeneDeletion = false;
        config.runDrugDeletion       = false;
        config.runFVA                = true;
        config.runFVASimilarity      = true;

end

%% ============================================================
% CREATE COMPARISON INFORMATION
% =============================================================

if numberOfTestModels == 1

    config.comparisonPairs = [];
    config.comparisonNames = {};

else

    config.comparisonPairs = [1, 2];

    config.comparisonNames = { ...
        [config.modelNames{1}, ...
         '_vs_', ...
         config.modelNames{2}]};

end

%% ============================================================
% CREATE TEST NAMES AND OUTPUT FOLDER
% =============================================================

config.datasetName = [ ...
    fullConfig.datasetName, ...
    '_', ...
    stageName];

config.outputPrefix = [ ...
    'TEST_', ...
    fullConfig.datasetName, ...
    '_', ...
    stageName];

config.outputFolder = fullfile( ...
    pwd, ...
    'results', ...
    config.outputPrefix);

config.saveIndividualGeneDeletionFiles = true;
config.saveIndividualDrugFiles         = true;
config.saveCompleteResults             = true;

%% ============================================================
% CHECK REQUIRED FUNCTIONS
% =============================================================

requiredFunctions = { ...
    'changeObjective', ...
    'optimizeCbModel', ...
    'singleGeneDeletion_rFASTCORMICS', ...
    'DrugDeletion', ...
    'fluxVariability', ...
    'FVA_similarity_Thomas'};

fprintf('\nChecking required MATLAB functions...\n');

for functionIndex = 1:numel(requiredFunctions)

    functionName = requiredFunctions{functionIndex};
    functionLocation = which(functionName);

    if isempty(functionLocation)

        error( ...
            'Required function not found: %s', ...
            functionName);

    end

    fprintf( ...
        '  Found: %-38s %s\n', ...
        functionName, ...
        functionLocation);

end

%% ============================================================
% CHECK SELECTED MODELS
% =============================================================

fprintf('\nChecking selected model structures...\n');

requiredModelFields = { ...
    'S', ...
    'rxns', ...
    'mets', ...
    'lb', ...
    'ub', ...
    'c'};

for modelIndex = 1:numel(config.models)

    model = config.models{modelIndex};
    modelName = config.modelNames{modelIndex};

    fprintf('\nModel %d: %s\n', ...
        modelIndex, modelName);

    for fieldIndex = 1:numel(requiredModelFields)

        fieldName = requiredModelFields{fieldIndex};

        if ~isfield(model, fieldName)

            error( ...
                'Model %s is missing field %s.', ...
                modelName, ...
                fieldName);

        end

    end

    numberOfReactions = numel(model.rxns);
    numberOfMetabolites = numel(model.mets);

    if size(model.S, 1) ~= numberOfMetabolites

        error( ...
            ['The number of rows in model.S does not match ' ...
             'model.mets for model %s.'], ...
            modelName);

    end

    if size(model.S, 2) ~= numberOfReactions

        error( ...
            ['The number of columns in model.S does not match ' ...
             'model.rxns for model %s.'], ...
            modelName);

    end

    if numel(model.lb) ~= numberOfReactions

        error( ...
            'model.lb size is incorrect for model %s.', ...
            modelName);

    end

    if numel(model.ub) ~= numberOfReactions

        error( ...
            'model.ub size is incorrect for model %s.', ...
            modelName);

    end

    if numel(model.c) ~= numberOfReactions

        error( ...
            'model.c size is incorrect for model %s.', ...
            modelName);

    end

    fprintf('  Reactions:   %d\n', numberOfReactions);
    fprintf('  Metabolites: %d\n', numberOfMetabolites);

    %% Check ATP objective

    atpFound = any(strcmp( ...
        model.rxns, ...
        config.atpObjective));

    if atpFound

        fprintf('  ATP objective found: %s\n', ...
            config.atpObjective);

    else

        warning( ...
            'ATP objective not found in model %s: %s', ...
            modelName, ...
            config.atpObjective);

    end

    %% Check biomass objective

    biomassFound = any(strcmp( ...
        model.rxns, ...
        config.biomassObjective));

    if biomassFound

        fprintf('  Biomass objective found: %s\n', ...
            config.biomassObjective);

    else

        error( ...
            'Biomass objective not found in model %s: %s', ...
            modelName, ...
            config.biomassObjective);

    end

    %% Check reaction mapping

    reactionFound = ismember( ...
        model.rxns, ...
        config.consistentModel.rxns);

    fprintf( ...
        '  Reference reactions mapped: %d of %d\n', ...
        sum(reactionFound), ...
        numberOfReactions);

    if any(~reactionFound)

        warning( ...
            ['Model %s has %d reactions that are not present ' ...
             'in consistent_model.'], ...
            modelName, ...
            sum(~reactionFound));

    end

    %% Test feasibility

    biomassModel = changeObjective( ...
        model, ...
        config.biomassObjective);

    solution = optimizeCbModel( ...
        biomassModel, ...
        'max');

    if isempty(solution) || ...
            ~isfield(solution, 'stat') || ...
            solution.stat ~= 1

        error( ...
            ['Model %s did not return an optimal ' ...
             'biomass solution.'], ...
            modelName);

    end

    fprintf( ...
        '  Biomass objective value: %.12g\n', ...
        solution.f);

end

%% ============================================================
% DISPLAY TEST SETTINGS
% =============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('Test configuration\n');
fprintf('============================================\n');

fprintf('Dataset: %s\n', fullConfig.datasetName);
fprintf('Stage:   %d\n', testStage);
fprintf('Name:    %s\n', stageName);

fprintf('Selected models:\n');

for modelIndex = 1:numel(config.modelNames)

    fprintf( ...
        '  %d. %s\n', ...
        modelIndex, ...
        config.modelNames{modelIndex});

end

fprintf('\nAnalyses:\n');
fprintf('  Single-gene deletion: %d\n', ...
    config.runSingleGeneDeletion);
fprintf('  Drug deletion:        %d\n', ...
    config.runDrugDeletion);
fprintf('  FVA:                  %d\n', ...
    config.runFVA);
fprintf('  FVA similarity:       %d\n', ...
    config.runFVASimilarity);

fprintf('\nOutput folder:\n');
fprintf('  %s\n', config.outputFolder);

%% ============================================================
% RUN THE UNIVERSAL ANALYSIS
% =============================================================

try

    results = universal_model_analysis(config);

catch analysisError

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('TEST FAILED\n');
    fprintf('============================================\n');

    fprintf('Identifier:\n');
    fprintf('  %s\n', analysisError.identifier);

    fprintf('Message:\n');
    fprintf('  %s\n', analysisError.message);

    fprintf('\nLocation:\n');

    for stackIndex = 1:numel(analysisError.stack)

        fprintf( ...
            '  %s | line %d\n', ...
            analysisError.stack(stackIndex).name, ...
            analysisError.stack(stackIndex).line);

    end

    rethrow(analysisError);

end

%% ============================================================
% VERIFY EXPECTED OUTPUTS
% =============================================================

fprintf('\nChecking analysis outputs...\n');

if config.runFVA

    if ~isfield(results, 'FVA')

        error( ...
            'The FVA stage ran but results.FVA was not created.');

    end

    for modelIndex = 1:numel(config.models)

        model = config.models{modelIndex};

        minFlux = ...
            results.FVA.minFluxModels{modelIndex};

        maxFlux = ...
            results.FVA.maxFluxModels{modelIndex};

        if numel(minFlux) ~= numel(model.rxns)

            error( ...
                ['Minimum FVA output size is incorrect ' ...
                 'for model %s.'], ...
                config.modelNames{modelIndex});

        end

        if numel(maxFlux) ~= numel(model.rxns)

            error( ...
                ['Maximum FVA output size is incorrect ' ...
                 'for model %s.'], ...
                config.modelNames{modelIndex});

        end

    end

    fprintf('  FVA output verified.\n');

end

if config.runDrugDeletion

    if ~isfield(results, 'drugDeletion')

        error( ...
            ['The drug-deletion stage ran but ' ...
             'results.drugDeletion was not created.']);

    end

    expectedDrugCount = numel( ...
        unique(config.GeneDrugRelations.DrugName));

    if size(results.drugDeletion.grRatio, 1) ~= ...
            expectedDrugCount

        error( ...
            'Drug-deletion output has an incorrect number of rows.');

    end

    if size(results.drugDeletion.grRatio, 2) ~= ...
            numel(config.models)

        error( ...
            'Drug-deletion output has an incorrect number of columns.');

    end

    fprintf('  Drug-deletion output verified.\n');

end

if config.runSingleGeneDeletion

    if ~isfield(results, 'singleGeneDeletion')

        error( ...
            ['The gene-deletion stage ran but ' ...
             'results.singleGeneDeletion was not created.']);

    end

    if numel(results.singleGeneDeletion) ~= ...
            numel(config.models)

        error( ...
            ['The number of gene-deletion result sets does not ' ...
             'match the number of models.']);

    end

    fprintf('  Single-gene deletion output verified.\n');

end

if config.runFVASimilarity

    if ~isfield(results, 'FVASimilarity')

        error( ...
            ['The FVA-similarity stage ran but ' ...
             'results.FVASimilarity was not created.']);

    end

    expectedNumberOfComparisons = ...
        size(config.comparisonPairs, 1);

    actualNumberOfComparisons = ...
        size( ...
            results.FVASimilarity.subsystemSimilarity, ...
            2);

    if actualNumberOfComparisons ~= ...
            expectedNumberOfComparisons

        error( ...
            ['The number of FVA similarity comparisons ' ...
             'is incorrect.']);

    end

    fprintf('  FVA-similarity output verified.\n');

end

%% ============================================================
% FINISHED
% =============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('TEST COMPLETED SUCCESSFULLY\n');
fprintf('============================================\n');

fprintf('Stage: %d - %s\n', ...
    testStage, stageName);

fprintf('Results saved in:\n');
fprintf('  %s\n', config.outputFolder);