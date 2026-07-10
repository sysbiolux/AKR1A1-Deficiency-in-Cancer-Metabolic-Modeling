clear;
close all;
clc;

%% MATLAB and COBRA information

fprintf('\n');
fprintf('============================================\n');
fprintf('Universal analysis test\n');
fprintf('MATLAB version: %s\n', version);
fprintf('============================================\n');

%% Check required functions

requiredFunctions = { ...
    'changeObjective', ...
    'singleGeneDeletion_rFASTCORMICS', ...
    'DrugDeletion', ...
    'fluxVariability', ...
    'FVA_similarity_Thomas'};

fprintf('\nChecking required functions...\n');

for functionIndex = 1:numel(requiredFunctions)

    functionName = requiredFunctions{functionIndex};

    functionLocation = which(functionName);

    if isempty(functionLocation)

        error( ...
            'Required function not found on the MATLAB path: %s', ...
            functionName);

    else

        fprintf('  Found %-40s %s\n', ...
            functionName, functionLocation);

    end

end

%% Select one dataset configuration
%
% Test the TCGA configuration first:
%
% fullConfig = config_TCGA_ONOFF();
%
% Or test the cell-line configuration:
%
% fullConfig = config_CellLines();

fullConfig = config_TCGA_ONOFF();

%% Select one model for the test
%
% TCGA examples:
%   1 = KIRC_ON
%   2 = KIRC_OFF
%   3 = KIRP_ON
%   4 = KIRP_OFF
%   5 = KICH_ON
%   6 = KICH_OFF
%   7 = LIHC_ON
%   8 = LIHC_OFF
%
% Cell-line examples:
%   1 = 769P_ctrl
%   2 = 769P_sc1
%   3 = 769P_sc2
%   4 = 769P_sc12
%   5 = Huh7_ctrl
%   6 = Huh7_sc1
%   7 = Huh7_sc2
%   8 = Huh7_sc12

testModelIndex = 1;

%% Verify the selected model index

numberOfAvailableModels = numel(fullConfig.models);

if testModelIndex < 1 || ...
        testModelIndex > numberOfAvailableModels

    error( ...
        'testModelIndex must be between 1 and %d.', ...
        numberOfAvailableModels);

end

%% Create a one-model configuration

config = fullConfig;

config.datasetName = [ ...
    fullConfig.datasetName, ...
    '_TEST_ONE_MODEL'];

config.models = ...
    fullConfig.models(testModelIndex);

config.modelNames = ...
    fullConfig.modelNames(testModelIndex);

if isfield(fullConfig, 'groupNames')

    config.groupNames = ...
        fullConfig.groupNames(testModelIndex);

end

if isfield(fullConfig, 'conditionNames')

    config.conditionNames = ...
        fullConfig.conditionNames(testModelIndex);

end

%% A single model has no between-model comparison

config.comparisonPairs = [];
config.comparisonNames = {};

%% Use a separate test output directory

config.outputFolder = fullfile( ...
    pwd, ...
    'Results_Test_One_Model');

config.outputPrefix = [ ...
    'TEST_', ...
    fullConfig.modelNames{testModelIndex}];

%% Select test stages
%
% Start with FVA only because it is usually the simplest test.
% Then activate drug deletion.
% Run single-gene deletion last because it can take much longer.

config.runSingleGeneDeletion = false;
config.runDrugDeletion       = false;
config.runFVA                = true;

% FVA similarity between models cannot be evaluated with one model.
config.runFVASimilarity      = false;

config.saveIndividualGeneDeletionFiles = true;
config.saveIndividualDrugFiles         = true;
config.saveCompleteResults             = true;

%% Display selected model

fprintf('\nSelected test model:\n');
fprintf('  Dataset: %s\n', fullConfig.datasetName);
fprintf('  Model index: %d\n', testModelIndex);
fprintf('  Model name: %s\n', config.modelNames{1});
fprintf('  Output folder: %s\n', config.outputFolder);

%% Validate the selected COBRA model

testModel = config.models{1};

requiredModelFields = { ...
    'S', ...
    'rxns', ...
    'mets', ...
    'lb', ...
    'ub', ...
    'c'};

fprintf('\nChecking COBRA model fields...\n');

for fieldIndex = 1:numel(requiredModelFields)

    fieldName = requiredModelFields{fieldIndex};

    if ~isfield(testModel, fieldName)

        error( ...
            'The selected model is missing the field: %s', ...
            fieldName);

    else

        fprintf('  Found model field: %s\n', fieldName);

    end

end

%% Check basic model dimensions

numberOfReactions = numel(testModel.rxns);
numberOfMetabolites = numel(testModel.mets);

if size(testModel.S, 1) ~= numberOfMetabolites

    error( ...
        ['The number of rows in model.S does not match ' ...
         'the number of metabolites.']);

end

if size(testModel.S, 2) ~= numberOfReactions

    error( ...
        ['The number of columns in model.S does not match ' ...
         'the number of reactions.']);

end

if numel(testModel.lb) ~= numberOfReactions

    error('The length of model.lb does not match model.rxns.');

end

if numel(testModel.ub) ~= numberOfReactions

    error('The length of model.ub does not match model.rxns.');

end

if numel(testModel.c) ~= numberOfReactions

    error('The length of model.c does not match model.rxns.');

end

fprintf('\nModel dimensions are consistent:\n');
fprintf('  Reactions: %d\n', numberOfReactions);
fprintf('  Metabolites: %d\n', numberOfMetabolites);

%% Check objectives

atpReactionIndex = find( ...
    strcmp(testModel.rxns, config.atpObjective));

biomassReactionIndex = find( ...
    strcmp(testModel.rxns, config.biomassObjective));

if isempty(atpReactionIndex)

    warning( ...
        'ATP objective reaction was not found: %s', ...
        config.atpObjective);

else

    fprintf('  ATP objective found: %s\n', ...
        config.atpObjective);

end

if isempty(biomassReactionIndex)

    warning( ...
        'Biomass objective reaction was not found: %s', ...
        config.biomassObjective);

else

    fprintf('  Biomass objective found: %s\n', ...
        config.biomassObjective);

end

%% Check mapping against the consistent model

[reactionFound, ~] = ismember( ...
    testModel.rxns, ...
    config.consistentModel.rxns);

numberMapped = sum(reactionFound);
numberNotMapped = sum(~reactionFound);

fprintf('\nReference-model mapping:\n');
fprintf('  Reactions mapped: %d\n', numberMapped);
fprintf('  Reactions not mapped: %d\n', numberNotMapped);

if numberNotMapped > 0

    warning( ...
        ['The selected model contains %d reactions that are not in ' ...
         'consistent_model.'], ...
        numberNotMapped);

end

%% Test model feasibility before the full analysis

fprintf('\nTesting model feasibility...\n');

biomassModel = changeObjective( ...
    testModel, ...
    config.biomassObjective);

solution = optimizeCbModel(biomassModel, 'max');

if isempty(solution) || ...
        ~isfield(solution, 'stat') || ...
        solution.stat ~= 1

    error( ...
        ['The selected model did not return an optimal biomass ' ...
         'solution. Check the solver and model constraints.']);

end

fprintf('  Biomass model is feasible.\n');
fprintf('  Objective value: %.12g\n', solution.f);

%% Run the universal function

fprintf('\nStarting the universal analysis test...\n');

try

    results = universal_model_analysis(config);

catch analysisError

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('UNIVERSAL ANALYSIS TEST FAILED\n');
    fprintf('============================================\n');

    fprintf('Error identifier:\n');
    fprintf('  %s\n', analysisError.identifier);

    fprintf('Error message:\n');
    fprintf('  %s\n', analysisError.message);

    fprintf('\nError location:\n');

    for stackIndex = 1:numel(analysisError.stack)

        fprintf( ...
            '  File: %s | Function: %s | Line: %d\n', ...
            analysisError.stack(stackIndex).file, ...
            analysisError.stack(stackIndex).name, ...
            analysisError.stack(stackIndex).line);

    end

    rethrow(analysisError);

end

%% Confirm expected results

if config.runFVA

    if ~isfield(results, 'FVA')
        error('The analysis finished but results.FVA was not created.');
    end

    if isempty(results.FVA.minFluxModels{1})
        error('The minimum FVA flux result is empty.');
    end

    if isempty(results.FVA.maxFluxModels{1})
        error('The maximum FVA flux result is empty.');
    end

    if numel(results.FVA.minFluxModels{1}) ~= numberOfReactions

        error( ...
            ['The number of minimum FVA values does not match ' ...
             'the number of reactions.']);

    end

    if numel(results.FVA.maxFluxModels{1}) ~= numberOfReactions

        error( ...
            ['The number of maximum FVA values does not match ' ...
             'the number of reactions.']);

    end

    fprintf('\nFVA output verified successfully.\n');

end

%% Final message

fprintf('\n');
fprintf('============================================\n');
fprintf('ONE-MODEL TEST COMPLETED SUCCESSFULLY\n');
fprintf('Model: %s\n', config.modelNames{1});
fprintf('============================================\n');