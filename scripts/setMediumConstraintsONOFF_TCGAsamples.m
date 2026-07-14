%% setMediumConstraintsONOFF_TCGAsamples.m
%
% Apply medium constraints to TCGA AKR1A1 ON/OFF models.
%
% Kidney cancers:
%   KICH, KIRC, KIRP -> RPMI
%
% Liver cancer:
%   LIHC -> DMEM
%
% MATLAB R2019b compatible.

clear;
close all;
clc;

%% ============================================================
% RELATIVE PATHS
%% ============================================================

modelDir = fullfile('results','models');
outputDir = fullfile('results','medium_constraints');
figureDir = fullfile(outputDir,'figures');

if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

if exist(figureDir,'dir') ~= 7
    mkdir(figureDir);
end

%% ============================================================
% LOAD TCGA MODELS
%% ============================================================

load(fullfile(modelDir,'driverModel_WithoutO2onoffLIHC.mat'), ...
    'model_LIHC_ON','model_LIHC_OFF');

load(fullfile(modelDir,'driverModel_WithoutO2onoffKIRC.mat'), ...
    'model_KIRC_ON','model_KIRC_OFF');

load(fullfile(modelDir,'driverModel_WithoutO2onoffKICH.mat'), ...
    'model_KICH_ON','model_KICH_OFF');

load(fullfile(modelDir,'driverModel_WithoutO2onoffKIRP.mat'), ...
    'model_KIRP_ON','model_KIRP_OFF');

%% ============================================================
% ORGANIZE MODELS
%% ============================================================

models_Kidney = {
    model_KIRC_ON
    model_KIRC_OFF
    model_KIRP_ON
    model_KIRP_OFF
    model_KICH_ON
    model_KICH_OFF
};

models_Liver = {
    model_LIHC_ON
    model_LIHC_OFF
};

modelNames_Kidney = {
    'model_KIRC_ON'
    'model_KIRC_OFF'
    'model_KIRP_ON'
    'model_KIRP_OFF'
    'model_KICH_ON'
    'model_KICH_OFF'
};

modelNames_Liver = {
    'model_LIHC_ON'
    'model_LIHC_OFF'
};

%% ============================================================
% RPMI MEDIUM
%% ============================================================

medium_composition_RPMI = {
    'EX_gly[e]'
    'EX_ala_L[e]'
    'EX_arg_L[e]'
    'EX_asn_L[e]'
    'EX_asp_L[e]'
    'EX_cys_L[e]'
    'EX_glu_L[e]'
    'EX_his_L[e]'
    'EX_4HPRO[e]'
    'EX_ile_L[e]'
    'EX_leu_L[e]'
    'EX_lys_L[e]'
    'EX_met_L[e]'
    'EX_phe_L[e]'
    'EX_pro_L[e]'
    'EX_ser_L[e]'
    'EX_thr_L[e]'
    'EX_trp_L[e]'
    'EX_tyr_L[e]'
    'EX_val_L[e]'
    'EX_btn[e]'
    'EX_chol[e]'
    'EX_pnto_R[e]'
    'EX_fol[e]'
    'EX_ncam[e]'
    'EX_pydxn[e]'
    'EX_ribflv[e]'
    'EX_thm[e]'
    'EX_inost[e]'
    'EX_ca2[e]'
    'EX_so4[e]'
    'EX_k[e]'
    'EX_hco3[e]'
    'EX_na1[e]'
    'EX_pi[e]'
    'EX_glc_D[e]'
    'EX_gln_L[e]'
    'EX_gthrd[e]'
};

met_Conc_mM_RPMI = [
    0.13333334
    2.0552995
    1.1494253
    0.37878788
    0.15037593
    0.20833333
    0.13605443
    0.09677419
    0.15267175
    0.3816794
    0.3816794
    0.21857923
    0.10067114
    0.09090909
    0.17391305
    0.2857143
    0.16806723
    0.024509804
    0.110497236
    0.17094018
    8.20E-04
    0.021428572
    0.000524
    0.002267574
    0.008196721
    0.004854369
    0.000532
    0.002967359
    0.19444445
    0.42372882
    0.40650406
    5.3333335
    23.809525
    103.44827
    5.633803
    11.111111
    2.0552995
    0.0032573289
];

%% ============================================================
% DMEM MEDIUM
%% ============================================================

medium_composition_DMEM = {
    'EX_gly[e]'
    'EX_ala_L[e]'
    'EX_arg_L[e]'
    'EX_cys_L[e]'
    'EX_his_L[e]'
    'EX_ile_L[e]'
    'EX_leu_L[e]'
    'EX_lys_L[e]'
    'EX_met_L[e]'
    'EX_phe_L[e]'
    'EX_ser_L[e]'
    'EX_thr_L[e]'
    'EX_trp_L[e]'
    'EX_tyr_L[e]'
    'EX_val_L[e]'
    'EX_chol[e]'
    'EX_pnto_R[e]'
    'EX_fol[e]'
    'EX_ncam[e]'
    'EX_pydxn[e]'
    'EX_ribflv[e]'
    'EX_thm[e]'
    'EX_inost[e]'
    'EX_ca2[e]'
    'EX_fe3[e]'
    'EX_so4[e]'
    'EX_k[e]'
    'EX_hco3[e]'
    'EX_na1[e]'
    'EX_pi[e]'
    'EX_glc_D[e]'
    'EX_pyr[e]'
    'EX_gln_L[e]'
};

met_Conc_mM_DMEM = [
    0.4
    3.9723501
    0.39810428
    0.20127796
    0.2
    0.8015267
    0.8015267
    0.7978142
    0.20134228
    0.4
    0.4
    0.79831934
    0.078431375
    0.39846742
    0.8034188
    0.028571429
    0.008385744
    0.009070295
    0.032786883
    0.019417476
    0.00106383
    0.011869436
    0.04
    1.8018018
    2.48E-04
    0.8139166
    5.3333335
    44.04762
    110.344826
    0.9057971
    25
    1
    3.9723501
];

%% Check medium dimensions

if numel(medium_composition_RPMI) ~= numel(met_Conc_mM_RPMI)
    error('RPMI reaction and concentration lists have different lengths.');
end

if numel(medium_composition_DMEM) ~= numel(met_Conc_mM_DMEM)
    error('DMEM reaction and concentration lists have different lengths.');
end

%% ============================================================
% APPLY RPMI TO KIDNEY MODELS
%% ============================================================

models_Kidney_RPMI = cell(size(models_Kidney));

for i = 1:numel(models_Kidney)

    model = models_Kidney{i};
    modelName = modelNames_Kidney{i};

    fprintf('\nApplying RPMI medium to %s\n',modelName);

    for ii = 1:numel(medium_composition_RPMI)

        rxn = medium_composition_RPMI{ii};
        conc = -met_Conc_mM_RPMI(ii);

        rxnIndex = find(strcmp(model.rxns,rxn));

        if ~isempty(rxnIndex)

            model.lb(rxnIndex) = conc;

        else

            fprintf('Reaction not found in %s: %s\n', ...
                modelName,rxn);
        end
    end

    models_Kidney_RPMI{i} = model;

    glucoseIndex = strcmp(model.rxns,'EX_glc_D[e]');

    if any(glucoseIndex)
        fprintf('Glucose lower bound: %g\n', ...
            model.lb(glucoseIndex));
    end
end

save(fullfile(outputDir,'models_Kidney_RPMI.mat'), ...
    'models_Kidney_RPMI');

%% ============================================================
% APPLY DMEM TO LIVER MODELS
%% ============================================================

models_Liver_DMEM = cell(size(models_Liver));

for i = 1:numel(models_Liver)

    model = models_Liver{i};
    modelName = modelNames_Liver{i};

    fprintf('\nApplying DMEM medium to %s\n',modelName);

    for ii = 1:numel(medium_composition_DMEM)

        rxn = medium_composition_DMEM{ii};
        conc = -met_Conc_mM_DMEM(ii);

        rxnIndex = find(strcmp(model.rxns,rxn));

        if ~isempty(rxnIndex)

            model.lb(rxnIndex) = conc;

        else

            fprintf('Reaction not found in %s: %s\n', ...
                modelName,rxn);
        end
    end

    models_Liver_DMEM{i} = model;

    glucoseIndex = strcmp(model.rxns,'EX_glc_D[e]');
    glutamineIndex = strcmp(model.rxns,'EX_gln_L[e]');

    if any(glucoseIndex)
        fprintf('Glucose lower bound: %g\n', ...
            model.lb(glucoseIndex));
    end

    if any(glutamineIndex)
        fprintf('Glutamine lower bound: %g\n', ...
            model.lb(glutamineIndex));
    end
end

save(fullfile(outputDir,'models_Liver_DMEM.mat'), ...
    'models_Liver_DMEM');

%% ============================================================
% INITIALIZE RESULTS
%% ============================================================

ResultsLiver = struct();
TFluxLiver = struct();

ResultsKidney = struct();
TFluxKidney = struct();

models_Liver_DMEMbiomass = cell(size(models_Liver_DMEM));
models_Kidney_RPMIbiomass = cell(size(models_Kidney_RPMI));

%% ============================================================
% LIVER FBA AND FVA
%% ============================================================

changeCobraSolver('glpk');

for i = 1:numel(models_Liver_DMEM)

    model = models_Liver_DMEM{i};
    modelName = modelNames_Liver{i};

    fprintf('\nRunning liver FBA/FVA: %s\n',modelName);

    model.lb(strcmp(model.rxns,'EX_HC02161[e]')) = 0;
    model.lb(strcmp(model.rxns,'EX_peplys[e]')) = 0;
    model.lb(strcmp(model.rxns,'EX_o2s[e]')) = 0;

    model.lb(strcmp(model.rxns,'EX_o2[e]')) = -50;
    model.lb(strcmp(model.rxns,'EX_co2[e]')) = -50;

    model = changeObjective(model,'biomass_reaction');

    FBALiver = optimizeCbModel(model,'max','zero');

    fprintf('Biomass objective value: %g\n',FBALiver.f);

    [minFluxLiver,maxFluxLiver] = fluxVariability(model,90);

    TFlux = table( ...
        model.rxns, ...
        FBALiver.x, ...
        minFluxLiver, ...
        maxFluxLiver, ...
        'VariableNames', ...
        {'Reaction','FBAFlux','MinFlux','MaxFlux'});

    TFluxLiver.(modelName) = TFlux;

    exchangeMask = contains(model.rxns,'EX_');

    TFluxEx = TFlux(exchangeMask,:);
    TFluxEx = sortrows(TFluxEx,'FBAFlux','ascend');

    ResultsLiver.(modelName) = TFluxEx;
    models_Liver_DMEMbiomass{i} = model;
end

save(fullfile(outputDir,'models_Liver_DMEMbiomass.mat'), ...
    'models_Liver_DMEMbiomass');

save(fullfile(outputDir,'ResultsLiver.mat'), ...
    'ResultsLiver','TFluxLiver');

%% ============================================================
% KIDNEY FBA AND FVA
%% ============================================================

for i = 1:numel(models_Kidney_RPMI)

    model = models_Kidney_RPMI{i};
    modelName = modelNames_Kidney{i};

    fprintf('\nRunning kidney FBA/FVA: %s\n',modelName);

    model.lb(strcmp(model.rxns,'EX_HC02161[e]')) = 0;
    model.lb(strcmp(model.rxns,'EX_peplys[e]')) = 0;
    model.lb(strcmp(model.rxns,'EX_o2s[e]')) = 0;

    model.lb(strcmp(model.rxns,'EX_o2[e]')) = -22.22;
    model.lb(strcmp(model.rxns,'EX_co2[e]')) = -22.22;

    model = changeObjective(model,'biomass_reaction');

    FBAKidney = optimizeCbModel(model,'max','zero');

    fprintf('Biomass objective value: %g\n',FBAKidney.f);

    [minFluxKidney,maxFluxKidney] = fluxVariability(model,90);

    TFlux = table( ...
        model.rxns, ...
        FBAKidney.x, ...
        minFluxKidney, ...
        maxFluxKidney, ...
        'VariableNames', ...
        {'Reaction','FBAFlux','MinFlux','MaxFlux'});

    TFluxKidney.(modelName) = TFlux;

    exchangeMask = contains(model.rxns,'EX_');

    TFluxEx = TFlux(exchangeMask,:);
    TFluxEx = sortrows(TFluxEx,'FBAFlux','ascend');

    ResultsKidney.(modelName) = TFluxEx;
    models_Kidney_RPMIbiomass{i} = model;
end

save(fullfile(outputDir,'models_Kidney_RPMIbiomass.mat'), ...
    'models_Kidney_RPMIbiomass');

save(fullfile(outputDir,'ResultsKidney.mat'), ...
    'ResultsKidney','TFluxKidney');

%% ============================================================
% SAVE COMPLETE WORKSPACE
%% ============================================================

save(fullfile(outputDir,'setMediumConstraintsONOFF_TCGAsamples.mat'));

%% ============================================================
% SELECT EXCHANGE REACTIONS FOR PLOTS
%% ============================================================

desiredReactions = {
    'EX_hco3[e]'
    'EX_glc_D[e]'
    'EX_o2[e]'
    'EX_pyr[e]'
    'EX_gln_L[e]'
    'EX_co2[e]'
    'EX_lac_L[e]'
    'EX_h2o[e]'
    'EX_ac[e]'
};

%% ============================================================
% KIDNEY EXCHANGE-FLUX PLOTS
%% ============================================================

modelsFluxEX = struct();

for i = 1:numel(models_Kidney_RPMIbiomass)

    modelName = modelNames_Kidney{i};
    TFluxEx = ResultsKidney.(modelName);

    selectedRows = ...
        ismember(TFluxEx.Reaction,desiredReactions) & ...
        TFluxEx.FBAFlux ~= 0;

    nonZeroTFluxEx = TFluxEx(selectedRows,:);

    if isempty(nonZeroTFluxEx)

        fprintf('No selected non-zero exchange fluxes for %s\n', ...
            modelName);

        continue
    end

    modelsFluxEX(i).flux = nonZeroTFluxEx;
    modelsFluxEX(i).modelname = modelName;

    figure;

    bar(nonZeroTFluxEx.FBAFlux);

    title(['Exchange Reaction Fluxes for ' modelName], ...
        'Interpreter','none');

    xlabel('Exchange reaction');
    ylabel('Flux');

    ax = gca;
    ax.XTick = 1:height(nonZeroTFluxEx);
    ax.XTickLabel = nonZeroTFluxEx.Reaction;

    xtickangle(90);

    saveas(gcf,fullfile(figureDir, ...
        [modelName '_EX_histogram.png']));

    saveas(gcf,fullfile(figureDir, ...
        [modelName '_EX_histogram.fig']));

    close(gcf);
end

%% ============================================================
% LIVER EXCHANGE-FLUX PLOTS
%% ============================================================

for i = 1:numel(models_Liver_DMEMbiomass)

    modelName = modelNames_Liver{i};
    TFluxEx = ResultsLiver.(modelName);

    selectedRows = ...
        ismember(TFluxEx.Reaction,desiredReactions) & ...
        TFluxEx.FBAFlux ~= 0;

    nonZeroTFluxEx = TFluxEx(selectedRows,:);

    if isempty(nonZeroTFluxEx)

        fprintf('No selected non-zero exchange fluxes for %s\n', ...
            modelName);

        continue
    end

    figure;

    bar(nonZeroTFluxEx.FBAFlux);

    title(['Non-Zero Exchange Reaction Fluxes for ' modelName], ...
        'Interpreter','none');

    xlabel('Exchange reaction');
    ylabel('Flux');

    ax = gca;
    ax.XTick = 1:height(nonZeroTFluxEx);
    ax.XTickLabel = nonZeroTFluxEx.Reaction;

    xtickangle(90);

    saveas(gcf,fullfile(figureDir, ...
        [modelName '_nonzero_histogram.png']));

    saveas(gcf,fullfile(figureDir, ...
        [modelName '_nonzero_histogram.fig']));

    close(gcf);
end

disp('setMediumConstraintsONOFF_TCGAsamples finished');