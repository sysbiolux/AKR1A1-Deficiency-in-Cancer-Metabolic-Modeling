function results = universal_model_analysis(config)
%UNIVERSAL_MODEL_ANALYSIS Run the same analysis for any model dataset.
%
% MATLAB version:
%   R2019b compatible
%
% This function does not change the scientific algorithms. It applies the
% existing analysis functions to every model supplied in config.models.
%
% Required input fields:
%
%   config.datasetName
%   config.models
%   config.modelNames
%   config.consistentModel
%   config.GeneDrugRelations
%   config.atpObjective
%   config.biomassObjective
%   config.deletionMethod
%   config.drugMethod
%   config.outputFolder
%   config.outputPrefix
%
% Optional execution fields:
%
%   config.runSingleGeneDeletion
%   config.runDrugDeletion
%   config.runFVA
%   config.runFVASimilarity
%
% Output:
%
%   results.singleGeneDeletion
%   results.drugDeletion
%   results.FVA
%   results.FVASimilarity

%% Validate the configuration

validateConfiguration(config);

numberOfModels = numel(config.models);

%% Add default execution settings when they were not supplied

if ~isfield(config, 'runSingleGeneDeletion')
    config.runSingleGeneDeletion = true;
end

if ~isfield(config, 'runDrugDeletion')
    config.runDrugDeletion = true;
end

if ~isfield(config, 'runFVA')
    config.runFVA = true;
end

if ~isfield(config, 'runFVASimilarity')
    config.runFVASimilarity = true;
end

if ~isfield(config, 'saveIndividualGeneDeletionFiles')
    config.saveIndividualGeneDeletionFiles = true;
end

if ~isfield(config, 'saveIndividualDrugFiles')
    config.saveIndividualDrugFiles = true;
end

if ~isfield(config, 'saveCompleteResults')
    config.saveCompleteResults = true;
end

%% Create the output folder

if ~exist(config.outputFolder, 'dir')
    mkdir(config.outputFolder);
end

%% Initialize the result structure

results = struct();

results.datasetName = config.datasetName;
results.modelNames  = config.modelNames;

if isfield(config, 'groupNames')
    results.groupNames = config.groupNames;
end

if isfield(config, 'conditionNames')
    results.conditionNames = config.conditionNames;
end

if isfield(config, 'comparisonPairs')
    results.comparisonPairs = config.comparisonPairs;
end

if isfield(config, 'comparisonNames')
    results.comparisonNames = config.comparisonNames;
end

%% Display dataset information

fprintf('\n');
fprintf('============================================\n');
fprintf('Dataset: %s\n', config.datasetName);
fprintf('Number of models: %d\n', numberOfModels);
fprintf('Output folder: %s\n', config.outputFolder);
fprintf('============================================\n');

%% Run single-gene deletion

if config.runSingleGeneDeletion

    fprintf('\nRunning single-gene deletion...\n');

    singleGeneDeletionResults = cell(numberOfModels, 1);

    for modelIndex = 1:numberOfModels

        modelName = config.modelNames{modelIndex};
        model = config.models{modelIndex};

        fprintf( ...
            '  Model %d of %d: %s\n', ...
            modelIndex, ...
            numberOfModels, ...
            modelName);

        %% ATP objective

        modelATP = changeObjective( ...
            model, ...
            config.atpObjective);

        [ ...
            grRatio_ATP, ...
            grRateKO_ATP, ...
            grRateWT_ATP, ...
            hasEffect_ATP, ...
            delRxns_ATP, ...
            fluxSolution_ATP, ...
            geneList_ATP ...
        ] = singleGeneDeletion_rFASTCORMICS( ...
            modelATP, ...
            config.deletionMethod, ...
            [], ...
            0, ...
            1);

        grRateWT_ATP = ...
            ones(numel(grRatio_ATP), 1) * grRateWT_ATP;

        %% Biomass objective

        modelBiomass = changeObjective( ...
            model, ...
            config.biomassObjective);

        [ ...
            grRatio_biomass, ...
            grRateKO_biomass, ...
            grRateWT_biomass, ...
            hasEffect_biomass, ...
            delRxns_biomass, ...
            fluxSolution_biomass, ...
            geneList_biomass ...
        ] = singleGeneDeletion_rFASTCORMICS( ...
            modelBiomass, ...
            config.deletionMethod, ...
            [], ...
            0, ...
            1);

        grRateWT_biomass = ...
            ones(numel(grRatio_biomass), 1) * grRateWT_biomass;

        %% Confirm that both analyses returned the same genes

        if numel(geneList_ATP) ~= numel(geneList_biomass)

            error( ...
                ['ATP and biomass deletion analyses returned different ' ...
                 'numbers of genes for model %s.'], ...
                modelName);

        end

        if ~isequal(geneList_ATP, geneList_biomass)

            warning( ...
                ['ATP and biomass gene-list order is different for ' ...
                 'model %s. The biomass gene list will be used.'], ...
                modelName);

        end

        %% Create the output table

        TEssential = table( ...
            geneList_biomass, ...
            grRatio_biomass, ...
            grRateKO_biomass, ...
            grRateWT_biomass, ...
            grRatio_ATP, ...
            grRateWT_ATP, ...
            grRateKO_ATP, ...
            'VariableNames', { ...
                'geneList', ...
                'grRatio_biomass', ...
                'grRateKO_biomass', ...
                'grRateWT_biomass', ...
                'grRatio_ATP', ...
                'grRateWT_ATP', ...
                'grRateKO_ATP'});

        %% Store complete details

        modelResult = struct();

        modelResult.modelName = modelName;
        modelResult.TEssential = TEssential;

        modelResult.hasEffect_ATP       = hasEffect_ATP;
        modelResult.delRxns_ATP         = delRxns_ATP;
        modelResult.fluxSolution_ATP    = fluxSolution_ATP;

        modelResult.hasEffect_biomass    = hasEffect_biomass;
        modelResult.delRxns_biomass      = delRxns_biomass;
        modelResult.fluxSolution_biomass = fluxSolution_biomass;

        singleGeneDeletionResults{modelIndex} = modelResult;

        %% Save one file for each model

        if config.saveIndividualGeneDeletionFiles

            safeModelName = makeSafeName(modelName);

            outputFile = fullfile( ...
                config.outputFolder, ...
                [config.outputPrefix, ...
                 '_essential_', ...
                 safeModelName, ...
                 '.mat']);

            save(outputFile, 'TEssential');

        end

    end

    results.singleGeneDeletion = singleGeneDeletionResults;

end

%% Run drug deletion

if config.runDrugDeletion

    fprintf('\nRunning drug deletion...\n');

    DrugList = unique( ...
        config.GeneDrugRelations.DrugName);

    numberOfDrugs = numel(DrugList);

    Drug_grRatio_biomass = ...
        zeros(numberOfDrugs, numberOfModels);

    Drug_grRateKO_biomass = ...
        zeros(numberOfDrugs, numberOfModels);

    Drug_grRateWT_biomass = ...
        zeros(numberOfDrugs, numberOfModels);

    lethalDrugs = cell(numberOfModels, 1);

    for modelIndex = 1:numberOfModels

        modelName = config.modelNames{modelIndex};
        model = config.models{modelIndex};

        fprintf( ...
            '  Model %d of %d: %s\n', ...
            modelIndex, ...
            numberOfModels, ...
            modelName);

        model = changeObjective( ...
            model, ...
            config.biomassObjective);

        [grRatio, grRateKO, grRateWT] = DrugDeletion( ...
            model, ...
            config.drugMethod, ...
            DrugList);

        %% Check output sizes before storing

        if numel(grRatio) ~= numberOfDrugs
            error( ...
                ['DrugDeletion returned an unexpected grRatio size ' ...
                 'for model %s.'], ...
                modelName);
        end

        if numel(grRateKO) ~= numberOfDrugs
            error( ...
                ['DrugDeletion returned an unexpected grRateKO size ' ...
                 'for model %s.'], ...
                modelName);
        end

        %% Store outputs

        Drug_grRatio_biomass(:, modelIndex) = grRatio;
        Drug_grRateKO_biomass(:, modelIndex) = grRateKO;

        if isscalar(grRateWT)

            Drug_grRateWT_biomass(:, modelIndex) = ...
                ones(numberOfDrugs, 1) * grRateWT;

        elseif numel(grRateWT) == numberOfDrugs

            Drug_grRateWT_biomass(:, modelIndex) = grRateWT;

        else

            error( ...
                ['DrugDeletion returned an unexpected grRateWT size ' ...
                 'for model %s.'], ...
                modelName);

        end

        lethalDrugs{modelIndex} = ...
            DrugList(grRatio == 0);

        %% Save lethal drugs for each model

        if config.saveIndividualDrugFiles

            safeModelName = makeSafeName(modelName);

            modelLethalDrugs = lethalDrugs{modelIndex};

            outputFile = fullfile( ...
                config.outputFolder, ...
                [config.outputPrefix, ...
                 '_lethalDrugs_', ...
                 safeModelName, ...
                 '.mat']);

            save(outputFile, 'modelLethalDrugs');

        end

    end

    %% Store drug-deletion outputs

    drugDeletionResult = struct();

    drugDeletionResult.DrugList = DrugList;

    drugDeletionResult.grRatio = ...
        Drug_grRatio_biomass;

    drugDeletionResult.grRateKO = ...
        Drug_grRateKO_biomass;

    drugDeletionResult.grRateWT = ...
        Drug_grRateWT_biomass;

    drugDeletionResult.lethalDrugs = ...
        lethalDrugs;

    results.drugDeletion = drugDeletionResult;

end

%% Run flux variability analysis

if config.runFVA

    fprintf('\nRunning flux variability analysis...\n');

    minFluxModels = cell(numberOfModels, 1);
    maxFluxModels = cell(numberOfModels, 1);
    fluxTables = cell(numberOfModels, 1);

    referenceRxns = config.consistentModel.rxns;
    numberOfReferenceRxns = numel(referenceRxns);

    minFluxReference = ...
        zeros(numberOfReferenceRxns, numberOfModels);

    maxFluxReference = ...
        zeros(numberOfReferenceRxns, numberOfModels);

    for modelIndex = 1:numberOfModels

        modelName = config.modelNames{modelIndex};
        model = config.models{modelIndex};

        fprintf( ...
            '  Model %d of %d: %s\n', ...
            modelIndex, ...
            numberOfModels, ...
            modelName);

        %% Preserve the current behavior:
        %  FVA uses the objective already present in each input model.

        [minFlux, maxFlux] = fluxVariability(model);

        minFluxModels{modelIndex} = minFlux;
        maxFluxModels{modelIndex} = maxFlux;

        %% Create the model-specific FVA table

        if isfield(model, 'rxnNames') && ...
                numel(model.rxnNames) == numel(model.rxns)

            TFlux = table( ...
                model.rxns, ...
                model.rxnNames, ...
                minFlux, ...
                maxFlux, ...
                'VariableNames', { ...
                    'rxns', ...
                    'rxnNames', ...
                    'minFlux', ...
                    'maxFlux'});

        else

            TFlux = table( ...
                model.rxns, ...
                minFlux, ...
                maxFlux, ...
                'VariableNames', { ...
                    'rxns', ...
                    'minFlux', ...
                    'maxFlux'});

        end

        fluxTables{modelIndex} = TFlux;

        %% Match model reactions to the consistent reference model

        [reactionFound, referenceLocation] = ...
            ismember(model.rxns, referenceRxns);

        if any(~reactionFound)

            warning( ...
                ['Model %s contains %d reactions that are not present ' ...
                 'in consistent_model. These reactions will not be ' ...
                 'included in the reference FVA matrices.'], ...
                modelName, ...
                sum(~reactionFound));

        end

        validModelRows = find(reactionFound);
        validReferenceRows = referenceLocation(reactionFound);

        minFluxReference(validReferenceRows, modelIndex) = ...
            minFlux(validModelRows);

        maxFluxReference(validReferenceRows, modelIndex) = ...
            maxFlux(validModelRows);

    end

    %% Store FVA outputs

    FVAResult = struct();

    FVAResult.minFluxModels = minFluxModels;
    FVAResult.maxFluxModels = maxFluxModels;
    FVAResult.fluxTables    = fluxTables;

    FVAResult.minFluxReference = minFluxReference;
    FVAResult.maxFluxReference = maxFluxReference;

    results.FVA = FVAResult;

end

%% Run FVA similarity

if config.runFVASimilarity

    if ~config.runFVA
        error('FVA must be enabled before FVA similarity can be calculated.');
    end

    fprintf('\nCalculating FVA similarity...\n');

    minFluxReference = results.FVA.minFluxReference;
    maxFluxReference = results.FVA.maxFluxReference;

    %% Calculate similarity for all model combinations

    allModelSimilarity = ...
        zeros(numberOfModels, numberOfModels);

    for firstModel = 1:numberOfModels

        for secondModel = 1:numberOfModels

            allModelSimilarity(firstModel, secondModel) = ...
                FVA_similarity_Thomas( ...
                    minFluxReference(:, firstModel), ...
                    maxFluxReference(:, firstModel), ...
                    minFluxReference(:, secondModel), ...
                    maxFluxReference(:, secondModel));

        end

    end

    %% Calculate only the comparisons specified in the configuration

    if isfield(config, 'comparisonPairs') && ...
            ~isempty(config.comparisonPairs)

        comparisonPairs = config.comparisonPairs;
        numberOfComparisons = size(comparisonPairs, 1);

    else

        comparisonPairs = createAllComparisonPairs(numberOfModels);
        numberOfComparisons = size(comparisonPairs, 1);

    end

    if isfield(config, 'comparisonNames') && ...
            numel(config.comparisonNames) == numberOfComparisons

        comparisonNames = config.comparisonNames(:);

    else

        comparisonNames = cell(numberOfComparisons, 1);

        for comparisonIndex = 1:numberOfComparisons

            firstModel = comparisonPairs(comparisonIndex, 1);
            secondModel = comparisonPairs(comparisonIndex, 2);

            comparisonNames{comparisonIndex} = sprintf( ...
                '%s_vs_%s', ...
                config.modelNames{firstModel}, ...
                config.modelNames{secondModel});

        end

    end

    %% Identify subsystems

    uniSys = unique( ...
        config.consistentModel.subSystems);

    numberOfSubsystems = numel(uniSys);

    subsystemSimilarity = ...
        zeros(numberOfSubsystems, numberOfComparisons);

    %% Calculate subsystem similarity

    for comparisonIndex = 1:numberOfComparisons

        firstModel = comparisonPairs(comparisonIndex, 1);
        secondModel = comparisonPairs(comparisonIndex, 2);

        fprintf( ...
            '  Comparison %d of %d: %s\n', ...
            comparisonIndex, ...
            numberOfComparisons, ...
            comparisonNames{comparisonIndex});

        for subsystemIndex = 1:numberOfSubsystems

            subsystemMatch = ismember( ...
                config.consistentModel.subSystems, ...
                uniSys(subsystemIndex));

            subsystemSimilarity( ...
                subsystemIndex, ...
                comparisonIndex) = ...
                FVA_similarity_Thomas( ...
                    minFluxReference(subsystemMatch, firstModel), ...
                    maxFluxReference(subsystemMatch, firstModel), ...
                    minFluxReference(subsystemMatch, secondModel), ...
                    maxFluxReference(subsystemMatch, secondModel));

        end

    end

    %% Create a readable table
    %
    % Each comparison gets its own table column.

    validComparisonNames = ...
        matlab.lang.makeValidName(comparisonNames);

    TFVASimilarity = array2table( ...
        subsystemSimilarity, ...
        'VariableNames', validComparisonNames);

    TFVASimilarity = addvars( ...
        TFVASimilarity, ...
        uniSys, ...
        'Before', ...
        1, ...
        'NewVariableNames', ...
        'Subsystem');

    %% Store outputs

    FVASimilarityResult = struct();

    FVASimilarityResult.allModelSimilarity = ...
        allModelSimilarity;

    FVASimilarityResult.comparisonPairs = ...
        comparisonPairs;

    FVASimilarityResult.comparisonNames = ...
        comparisonNames;

    FVASimilarityResult.uniSys = ...
        uniSys;

    FVASimilarityResult.subsystemSimilarity = ...
        subsystemSimilarity;

    FVASimilarityResult.table = ...
        TFVASimilarity;

    results.FVASimilarity = FVASimilarityResult;

    %% Save the readable similarity table

    similarityFile = fullfile( ...
        config.outputFolder, ...
        [config.outputPrefix, ...
         '_FVA_similarity.mat']);

    save( ...
        similarityFile, ...
        'TFVASimilarity', ...
        'allModelSimilarity', ...
        'comparisonPairs', ...
        'comparisonNames');

end

%% Save the complete result structure

if config.saveCompleteResults

    completeOutputFile = fullfile( ...
        config.outputFolder, ...
        [config.outputPrefix, ...
         '_complete_analysis.mat']);

    save( ...
        completeOutputFile, ...
        'results', ...
        'config', ...
        '-v7.3');

end

fprintf('\nAnalysis finished for dataset %s.\n', config.datasetName);

end


function validateConfiguration(config)
%VALIDATECONFIGURATION Validate required universal inputs.

requiredFields = { ...
    'datasetName', ...
    'models', ...
    'modelNames', ...
    'consistentModel', ...
    'GeneDrugRelations', ...
    'atpObjective', ...
    'biomassObjective', ...
    'deletionMethod', ...
    'drugMethod', ...
    'outputFolder', ...
    'outputPrefix'};

for fieldIndex = 1:numel(requiredFields)

    fieldName = requiredFields{fieldIndex};

    if ~isfield(config, fieldName)

        error( ...
            'The configuration is missing the required field: %s', ...
            fieldName);

    end

end

if ~iscell(config.models)
    error('config.models must be a cell array of COBRA models.');
end

if ~iscell(config.modelNames)
    error('config.modelNames must be a cell array of names.');
end

if numel(config.models) ~= numel(config.modelNames)

    error( ...
        ['The number of models (%d) does not match the number ' ...
         'of model names (%d).'], ...
        numel(config.models), ...
        numel(config.modelNames));

end

if isempty(config.models)
    error('No models were supplied in config.models.');
end

if ~isfield(config.consistentModel, 'rxns')
    error('config.consistentModel does not contain the field rxns.');
end

if ~isfield(config.consistentModel, 'subSystems')

    error( ...
        'config.consistentModel does not contain the field subSystems.');

end

if ~isfield(config.GeneDrugRelations, 'DrugName')

    error( ...
        ['config.GeneDrugRelations does not contain the ' ...
         'field DrugName.']);

end

if isfield(config, 'comparisonPairs') && ...
        ~isempty(config.comparisonPairs)

    comparisonPairs = config.comparisonPairs;
    numberOfModels = numel(config.models);

    if size(comparisonPairs, 2) ~= 2

        error( ...
            'config.comparisonPairs must have exactly two columns.');

    end

    if any(comparisonPairs(:) < 1) || ...
            any(comparisonPairs(:) > numberOfModels)

        error( ...
            ['config.comparisonPairs contains a model index outside ' ...
             'the valid range 1 to %d.'], ...
            numberOfModels);

    end

end

end


function safeName = makeSafeName(inputName)
%MAKESAFENAME Create a filename-safe model name.

safeName = regexprep(inputName, '[^a-zA-Z0-9_-]', '_');

end


function comparisonPairs = createAllComparisonPairs(numberOfModels)
%CREATEALLCOMPARISONPAIRS Return each unique model pair once.

numberOfPairs = ...
    numberOfModels * (numberOfModels - 1) / 2;

comparisonPairs = zeros(numberOfPairs, 2);

pairIndex = 0;

for firstModel = 1:numberOfModels

    for secondModel = 1:firstModel - 1

        pairIndex = pairIndex + 1;

        comparisonPairs(pairIndex, :) = [ ...
            firstModel, ...
            secondModel];

    end

end

end