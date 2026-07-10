%% setMediumConstraints_unified.m
%
% Applies media constraints to:
%
% TCGA
%   KICH, KIRC, KIRP -> RPMI
%   LIHC             -> DMEM
%
% CELL LINES
%   769-p            -> RPMI
%   Huh7             -> DMEM
%
% Models are loaded from:
%   cfg.modelDir
%
% Results are saved to:
%   cfg.modelsDir
%   cfg.analysisDir
%
% MATLAB R2019b compatible.

clear;
close all;
clc;

cfg = pipeline_config;

%% ============================================================
% CHECK OUTPUT DIRECTORIES
%% ============================================================

if exist(cfg.modelsDir, 'dir') ~= 7
    mkdir(cfg.modelsDir);
end

if exist(cfg.analysisDir, 'dir') ~= 7
    mkdir(cfg.analysisDir);
end

%% ============================================================
% INITIALIZE COBRA TOOLBOX
%% ============================================================

if exist('initCobraToolbox', 'file') == 2
    initCobraToolbox(false);
else
    warning('initCobraToolbox was not found on the MATLAB path.');
end

if exist('changeCobraSolver', 'file') ~= 2
    error('COBRA Toolbox function changeCobraSolver was not found.');
end

if exist('optimizeCbModel', 'file') ~= 2
    error('COBRA Toolbox function optimizeCbModel was not found.');
end

if exist('fluxVariability', 'file') ~= 2
    error('COBRA Toolbox function fluxVariability was not found.');
end

changeCobraSolver('glpk');

%% ============================================================
% MODEL DEFINITIONS
%% ============================================================

% Each row:
%   model variable name
%   MAT filename
%   biological group
%   medium

modelList = { ...
    'model_KICH_ON',   'model_KICH_ON.mat',   'TCGA_KIDNEY', 'RPMI'; ...
    'model_KICH_OFF',  'model_KICH_OFF.mat',  'TCGA_KIDNEY', 'RPMI'; ...
    'model_KIRC_ON',   'model_KIRC_ON.mat',   'TCGA_KIDNEY', 'RPMI'; ...
    'model_KIRC_OFF',  'model_KIRC_OFF.mat',  'TCGA_KIDNEY', 'RPMI'; ...
    'model_KIRP_ON',   'model_KIRP_ON.mat',   'TCGA_KIDNEY', 'RPMI'; ...
    'model_KIRP_OFF',  'model_KIRP_OFF.mat',  'TCGA_KIDNEY', 'RPMI'; ...
    'model_LIHC_ON',   'model_LIHC_ON.mat',   'TCGA_LIVER',  'DMEM'; ...
    'model_LIHC_OFF',  'model_LIHC_OFF.mat',  'TCGA_LIVER',  'DMEM'; ...
    'model_7_ctrl',    'model_7_ctrl.mat',    'CELL_769P',   'RPMI'; ...
    'model_7_sc1',     'model_7_sc1.mat',     'CELL_769P',   'RPMI'; ...
    'model_7_sc12',    'model_7_sc12.mat',    'CELL_769P',   'RPMI'; ...
    'model_7_sc2',     'model_7_sc2.mat',     'CELL_769P',   'RPMI'; ...
    'model_H_ctrl',    'model_H_ctrl.mat',    'CELL_HUH7',   'DMEM'; ...
    'model_H_sc1',     'model_H_sc1.mat',     'CELL_HUH7',   'DMEM'; ...
    'model_H_sc12',    'model_H_sc12.mat',    'CELL_HUH7',   'DMEM'; ...
    'model_H_sc2',     'model_H_sc2.mat',     'CELL_HUH7',   'DMEM'};

%% ============================================================
% RPMI MEDIUM
%% ============================================================

mediumRxnsRPMI = { ...
    'EX_gly[e]';
    'EX_ala_L[e]';
    'EX_arg_L[e]';
    'EX_asn_L[e]';
    'EX_asp_L[e]';
    'EX_cys_L[e]';
    'EX_glu_L[e]';
    'EX_his_L[e]';
    'EX_4HPRO[e]';
    'EX_ile_L[e]';
    'EX_leu_L[e]';
    'EX_lys_L[e]';
    'EX_met_L[e]';
    'EX_phe_L[e]';
    'EX_pro_L[e]';
    'EX_ser_L[e]';
    'EX_thr_L[e]';
    'EX_trp_L[e]';
    'EX_tyr_L[e]';
    'EX_val_L[e]';
    'EX_btn[e]';
    'EX_chol[e]';
    'EX_pnto_R[e]';
    'EX_fol[e]';
    'EX_ncam[e]';
    'EX_pydxn[e]';
    'EX_ribflv[e]';
    'EX_thm[e]';
    'EX_inost[e]';
    'EX_ca2[e]';
    'EX_so4[e]';
    'EX_k[e]';
    'EX_hco3[e]';
    'EX_na1[e]';
    'EX_pi[e]';
    'EX_glc_D[e]';
    'EX_gln_L[e]';
    'EX_gthrd[e]'};

mediumConcRPMI = [ ...
    0.13333334;
    2.0552995;
    1.1494253;
    0.37878788;
    0.15037593;
    0.20833333;
    0.13605443;
    0.09677419;
    0.15267175;
    0.3816794;
    0.3816794;
    0.21857923;
    0.10067114;
    0.09090909;
    0.17391305;
    0.2857143;
    0.16806723;
    0.024509804;
    0.110497236;
    0.17094018;
    8.20e-4;
    0.021428572;
    0.000524;
    0.002267574;
    0.008196721;
    0.004854369;
    0.000532;
    0.002967359;
    0.19444445;
    0.42372882;
    0.40650406;
    5.3333335;
    23.809525;
    103.44827;
    5.633803;
    11.111111;
    2.0552995;
    0.0032573289];

%% ============================================================
% DMEM MEDIUM
%% ============================================================

mediumRxnsDMEM = { ...
    'EX_gly[e]';
    'EX_ala_L[e]';
    'EX_arg_L[e]';
    'EX_cys_L[e]';
    'EX_his_L[e]';
    'EX_ile_L[e]';
    'EX_leu_L[e]';
    'EX_lys_L[e]';
    'EX_met_L[e]';
    'EX_phe_L[e]';
    'EX_ser_L[e]';
    'EX_thr_L[e]';
    'EX_trp_L[e]';
    'EX_tyr_L[e]';
    'EX_val_L[e]';
    'EX_chol[e]';
    'EX_pnto_R[e]';
    'EX_fol[e]';
    'EX_ncam[e]';
    'EX_pydxn[e]';
    'EX_ribflv[e]';
    'EX_thm[e]';
    'EX_inost[e]';
    'EX_ca2[e]';
    'EX_fe3[e]';
    'EX_so4[e]';
    'EX_k[e]';
    'EX_hco3[e]';
    'EX_na1[e]';
    'EX_pi[e]';
    'EX_glc_D[e]';
    'EX_pyr[e]';
    'EX_gln_L[e]'};

mediumConcDMEM = [ ...
    0.4;
    3.9723501;
    0.39810428;
    0.20127796;
    0.2;
    0.8015267;
    0.8015267;
    0.7978142;
    0.20134228;
    0.4;
    0.4;
    0.79831934;
    0.078431375;
    0.39846742;
    0.8034188;
    0.028571429;
    0.008385744;
    0.009070295;
    0.032786883;
    0.019417476;
    0.00106383;
    0.011869436;
    0.04;
    1.8018018;
    2.48e-4;
    0.8139166;
    5.3333335;
    44.04762;
    110.344826;
    0.9057971;
    25;
    1;
    3.9723501];

%% Validate medium definitions

if numel(mediumRxnsRPMI) ~= numel(mediumConcRPMI)
    error('RPMI reaction and concentration counts do not match.');
end

if numel(mediumRxnsDMEM) ~= numel(mediumConcDMEM)
    error('DMEM reaction and concentration counts do not match.');
end

%% ============================================================
% LOAD, CONSTRAIN AND ANALYZE EVERY MODEL
%% ============================================================

constrainedModels = struct();
FBAResults = struct();
FVAResults = struct();
ExchangeResults = struct();
missingReactions = struct();

for iModel = 1:size(modelList, 1)

    modelVariable = modelList{iModel, 1};
    modelFileName = modelList{iModel, 2};
    groupName = modelList{iModel, 3};
    mediumName = modelList{iModel, 4};

    modelFile = fullfile(cfg.modelDir, modelFileName);

    fprintf('\n========================================\n');
    fprintf('Processing: %s\n', modelVariable);
    fprintf('Group: %s\n', groupName);
    fprintf('Medium: %s\n', mediumName);
    fprintf('File: %s\n', modelFile);
    fprintf('========================================\n');

    if exist(modelFile, 'file') ~= 2
        warning('Model file not found. Skipping: %s', modelFile);
        continue;
    end

    %% Load model from individual MAT file

    loadedData = load(modelFile);

    if isfield(loadedData, modelVariable)

        model = loadedData.(modelVariable);

    elseif isfield(loadedData, 'ContextModel')

        model = loadedData.ContextModel;

    else

        availableVariables = fieldnames(loadedData);

        if numel(availableVariables) == 1
            model = loadedData.(availableVariables{1});
            warning( ...
                'Using variable "%s" from %s.', ...
                availableVariables{1}, ...
                modelFileName);
        else
            error( ...
                ['Cannot identify the model variable in:\n%s\n' ...
                 'Expected "%s" or "ContextModel".'], ...
                 modelFile, ...
                 modelVariable);
        end

    end

    if ~isstruct(model) || ~isfield(model, 'rxns')
        error('Loaded variable from %s is not a COBRA model.', modelFile);
    end

    %% Choose medium according to biological system

    if strcmpi(mediumName, 'RPMI')
        mediumRxns = mediumRxnsRPMI;
        mediumConc = mediumConcRPMI;
    elseif strcmpi(mediumName, 'DMEM')
        mediumRxns = mediumRxnsDMEM;
        mediumConc = mediumConcDMEM;
    else
        error('Unknown medium: %s', mediumName);
    end

    %% Apply medium lower bounds

    missingForModel = {};

    for iMedium = 1:numel(mediumRxns)

        rxn = mediumRxns{iMedium};
        lowerBound = -mediumConc(iMedium);

        rxnIndex = find(strcmp(model.rxns, rxn));

        if isempty(rxnIndex)
            missingForModel{end + 1, 1} = rxn;
        else
            model.lb(rxnIndex) = lowerBound;
        end

    end

    missingReactions.(modelVariable) = missingForModel;

    fprintf('Medium reactions applied: %d\n', ...
        numel(mediumRxns) - numel(missingForModel));

    fprintf('Medium reactions missing: %d\n', ...
        numel(missingForModel));

    %% Additional original constraints

    model = set_bound_if_present( ...
        model, 'EX_HC02161[e]', 0, 'l');

    model = set_bound_if_present( ...
        model, 'EX_peplys[e]', 0, 'l');

    model = set_bound_if_present( ...
        model, 'EX_o2s[e]', 0, 'l');

    if strcmpi(mediumName, 'RPMI')

        model = set_bound_if_present( ...
            model, 'EX_o2[e]', -22.22, 'l');

        model = set_bound_if_present( ...
            model, 'EX_co2[e]', -22.22, 'l');

    else

        model = set_bound_if_present( ...
            model, 'EX_o2[e]', -50, 'l');

        model = set_bound_if_present( ...
            model, 'EX_co2[e]', -50, 'l');

    end

    %% Set biomass objective

    if any(strcmp(model.rxns, 'biomass_reaction'))
        model = changeObjective(model, 'biomass_reaction');
    else
        warning( ...
            'biomass_reaction not found in %s.', ...
            modelVariable);
    end

    %% Run FBA

    FBASolution = optimizeCbModel(model, 'max', 'zero');

    FBAResults.(modelVariable) = FBASolution;

    if isempty(FBASolution.f)
        fprintf('FBA returned an empty objective value.\n');
    else
        fprintf('FBA objective: %.10g\n', FBASolution.f);
    end

    %% Run FVA

    [minFlux, maxFlux] = fluxVariability(model, 90);

    fluxTable = table( ...
        model.rxns, ...
        FBASolution.x, ...
        minFlux, ...
        maxFlux, ...
        'VariableNames', { ...
            'Reaction', ...
            'FBAFlux', ...
            'MinimumFlux', ...
            'MaximumFlux'});

    FVAResults.(modelVariable) = fluxTable;

    %% Extract exchange reactions

    exchangeMask = contains(model.rxns, 'EX_');

    exchangeTable = fluxTable(exchangeMask, :);

    exchangeTable = sortrows( ...
        exchangeTable, ...
        'FBAFlux', ...
        'ascend');

    ExchangeResults.(modelVariable) = exchangeTable;

    %% Keep constrained model

    constrainedModels.(modelVariable) = model;

    %% Save the individual constrained model

    constrainedModelFile = fullfile( ...
        cfg.modelsDir, ...
        [modelVariable '_constrained.mat']);

    constrainedModel = model;

    save( ...
        constrainedModelFile, ...
        'constrainedModel', ...
        'mediumName', ...
        'groupName');

    fprintf('Saved constrained model:\n%s\n', ...
        constrainedModelFile);

end

%% ============================================================
% CREATE GROUPED MODEL COLLECTIONS
%% ============================================================

models_Kidney_RPMI = {};
modelNames_Kidney = {};

models_Liver_DMEM = {};
modelNames_Liver = {};

models_769p_RPMI = {};
modelNames_769p = {};

models_Huh7_DMEM = {};
modelNames_Huh7 = {};

for iModel = 1:size(modelList, 1)

    modelVariable = modelList{iModel, 1};
    groupName = modelList{iModel, 3};

    if ~isfield(constrainedModels, modelVariable)
        continue;
    end

    switch groupName

        case 'TCGA_KIDNEY'

            models_Kidney_RPMI{end + 1} = ...
                constrainedModels.(modelVariable);

            modelNames_Kidney{end + 1} = ...
                modelVariable;

        case 'TCGA_LIVER'

            models_Liver_DMEM{end + 1} = ...
                constrainedModels.(modelVariable);

            modelNames_Liver{end + 1} = ...
                modelVariable;

        case 'CELL_769P'

            models_769p_RPMI{end + 1} = ...
                constrainedModels.(modelVariable);

            modelNames_769p{end + 1} = ...
                modelVariable;

        case 'CELL_HUH7'

            models_Huh7_DMEM{end + 1} = ...
                constrainedModels.(modelVariable);

            modelNames_Huh7{end + 1} = ...
                modelVariable;

    end

end

%% ============================================================
% SAVE GROUPED MODELS
%% ============================================================

save( ...
    fullfile(cfg.modelsDir, 'models_Kidney_RPMI.mat'), ...
    'models_Kidney_RPMI', ...
    'modelNames_Kidney');

save( ...
    fullfile(cfg.modelsDir, 'models_Liver_DMEM.mat'), ...
    'models_Liver_DMEM', ...
    'modelNames_Liver');

save( ...
    fullfile(cfg.modelsDir, 'models_769p_RPMI.mat'), ...
    'models_769p_RPMI', ...
    'modelNames_769p');

save( ...
    fullfile(cfg.modelsDir, 'models_Huh7_DMEM.mat'), ...
    'models_Huh7_DMEM', ...
    'modelNames_Huh7');

%% ============================================================
% SAVE ANALYSIS RESULTS
%% ============================================================

save( ...
    fullfile(cfg.analysisDir, ...
    'medium_constraint_FBA_FVA_results.mat'), ...
    'FBAResults', ...
    'FVAResults', ...
    'ExchangeResults', ...
    'missingReactions', ...
    'modelList', ...
    '-v7.3');

fprintf('\n========================================\n');
fprintf('Medium constraint analysis completed.\n');
fprintf('Constrained models: %s\n', cfg.modelsDir);
fprintf('Analysis results: %s\n', cfg.analysisDir);
fprintf('========================================\n');

%% ============================================================
% LOCAL FUNCTION
%% ============================================================

function model = set_bound_if_present( ...
    model, reactionID, boundValue, boundType)

reactionIndex = find(strcmp(model.rxns, reactionID));

if isempty(reactionIndex)
    return;
end

if strcmpi(boundType, 'l')
    model.lb(reactionIndex) = boundValue;
elseif strcmpi(boundType, 'u')
    model.ub(reactionIndex) = boundValue;
else
    error('boundType must be "l" or "u".');
end

end