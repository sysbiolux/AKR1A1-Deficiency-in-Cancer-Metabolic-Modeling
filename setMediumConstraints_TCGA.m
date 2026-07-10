%%
Initialize the COBRA Toolbox
initCobraToolbox(false); % Set true if you want to update the toolbox

load('driverModel_WithoutO2onoffLIHC','model_LIHC_ON', 'model_LIHC_OFF');                       
load('driverModel_WithoutO2onoffKIRC','model_KIRC_ON', 'model_KIRC_OFF');  
load('driverModel_WithoutO2onoffKICH','model_KICH_ON', 'model_KICH_OFF');  
load('driverModel_WithoutO2onoffKIRP','model_KIRP_ON', 'model_KIRP_OFF'); 

models_TCGA = {model_KIRC_ON, model_KIRC_OFF, model_KIRP_ON, model_KIRP_OFF, model_KICH_ON,...
               model_KICH_OFF, model_LIHC_ON, model_LIHC_OFF};

modelNames_TCGA = {'model_KIRC_ON', 'model_KIRC_OFF', 'model_KIRP_ON', 'model_KIRP_OFF', 'model_KICH_ON',...
                   'model_KICH_OFF', 'model_LIHC_ON', 'model_LIHC_OFF'};

model_constrained=struct();
%%
% Medium for tumor microenvironment
medium_composition_tumor = {'EX_glc__D_e'; 'EX_gln__L_e'; 'EX_lys__L_e'; 'EX_leu__L_e'; 'EX_met__L_e'; ...
                            'EX_phe__L_e'; 'EX_thr__L_e'; 'EX_trp__L_e'; 'EX_val__L_e'; 'EX_ala__L_e'; ...
                            'EX_ser__L_e'; 'EX_asp__L_e'; 'EX_glu__L_e'; 'EX_pro__L_e'; 'EX_palm_e'; ...
                            'EX_ole_e'; 'EX_stear_e'; 'EX_lac__L_e'; 'EX_o2_e'; 'EX_thm_e'; ...
                            'EX_ribflv_e'; 'EX_nac_e'};
                        
met_Conc_mM_tumor = [7.5; 1.3; 0.1; 0.15; 0.03; ...
                     0.06; 0.1; 0.05; 0.2; 0.35; ...
                     0.15; 0.02; 0.06; 0.2; 0.3; ...
                     0.2; 0.08; 10; 1; 0.001; ...
                     0.001; 0.001];


%% Loop over the medium compositions for TCGA models
for ii = 1:numel(medium_composition_TCGAmodels)
    rxn = medium_composition_tumor{ii};
    conc = -met_Conc_mM_tumor(ii);
    
    for i = 1:numel(models_TCGA)
        model = models_TCGA{i};
        modelName = modelNames_TCGA{i};
        %modelSampling = samplingResults(i).modelSampling;
        
        % Find the index of the reaction corresponding to the metabolite in the model
        rxnIndex = find(strcmp(model.rxns, rxn));
        if ~isempty(rxnIndex)
            model.lb(rxnIndex) =  conc;
            disp(model.lb(ismember(model.rxns,'EX_glc_D[e]')));
        else
            disp(['Reaction not found in ' modelName ': ' rxn]);
        end
        
        models_TCGA{i}= model;% save here every new model transformed with the addition of the media
        
    end
end

save models_TCGA models_TCGA

%%
%% 7690 7 
for i=1:numel(models_TCGA)
changeCobraSolver('glpk')
model=models_Kidney_TCGAmodels{i};
model.lb(ismember(model.rxns,'EX_HC02161[e]'))=0;
model.lb(ismember(model.rxns,'EX_peplys[e]'))=0;
%model.lb(ismember(model.rxns,'EX_hdl_hs[e]'))=0;
%model.lb(ismember(model.rxns,'EX_ldl_hs[e]'))=0;
model.lb(ismember(model.rxns,'EX_o2s[e]'))=0;
model.lb(ismember(model.rxns,'EX_o2[e]'))=-22.22; %Glc x2
model.lb(ismember(model.rxns,'EX_co2[e]'))=-22.22; %Glc x2
changeObjective(model,'biomass_reaction');
FBAKidney=optimizeCbModel(model,'max','zero');
FBAKidney.f
[minFluxKidney, maxFluxKidney] = fluxVariability(model, 90);
TFlux=table(model.rxns,FBAKidney.x,minFluxKidney,maxFluxKidney);
TFluxKidney.(modelNames_Kidney{i})=TFlux;
TFluxEx=TFlux(contains(model.rxns, 'EX_'),:);
TFluxEx = sortrows(TFluxEx,'Var2','ascend');
ResultsKidney.(modelNames_Kidney{i})=TFluxEx;
models_Kidney_TCGAbiomass{i}=model;
end

save models_TCGAbiomass models_TCGAbiomass

%%
% Load your GEM model
model = readCbModel('your_cancer_GEM.mat');

% Define media composition constraints
% Glucose
model = changeRxnBounds(model, 'EX_glc__D_e', -7.5, 'l'); % Midpoint of the 5-10 mM range

% L-Glutamine
model = changeRxnBounds(model, 'EX_gln__L_e', -1.3, 'l'); % Midpoint of the 0.6-2 mM range

% Essential Amino Acids
model = changeRxnBounds(model, 'EX_lys__L_e', -0.1, 'l'); % L-Lysine
model = changeRxnBounds(model, 'EX_leu__L_e', -0.15, 'l'); % L-Leucine
model = changeRxnBounds(model, 'EX_met__L_e', -0.03, 'l'); % L-Methionine
model = changeRxnBounds(model, 'EX_phe__L_e', -0.06, 'l'); % L-Phenylalanine
model = changeRxnBounds(model, 'EX_thr__L_e', -0.1, 'l'); % L-Threonine
model = changeRxnBounds(model, 'EX_trp__L_e', -0.05, 'l'); % L-Tryptophan
model = changeRxnBounds(model, 'EX_val__L_e', -0.2, 'l'); % L-Valine

% Non-Essential Amino Acids
model = changeRxnBounds(model, 'EX_ala__L_e', -0.35, 'l'); % L-Alanine
model = changeRxnBounds(model, 'EX_ser__L_e', -0.15, 'l'); % L-Serine
model = changeRxnBounds(model, 'EX_asp__L_e', -0.02, 'l'); % L-Aspartate
model = changeRxnBounds(model, 'EX_glu__L_e', -0.06, 'l'); % L-Glutamate
model = changeRxnBounds(model, 'EX_pro__L_e', -0.2, 'l'); % L-Proline

% Fatty Acids
model = changeRxnBounds(model, 'EX_palm_e', -0.3, 'l'); % Palmitate
model = changeRxnBounds(model, 'EX_ole_e', -0.2, 'l'); % Oleate
model = changeRxnBounds(model, 'EX_stear_e', -0.08, 'l'); % Stearate

% Lactate (allow export)
model = changeRxnBounds(model, 'EX_lac__L_e', 10, 'u'); % Lactate can be high due to Warburg effect

% Oxygen (hypoxia condition)
model = changeRxnBounds(model, 'EX_o2_e', -1, 'l'); % Adjust based on desired hypoxia level

% Vitamins and Minerals (example, adjust as necessary)
model = changeRxnBounds(model, 'EX_thm_e', -0.001, 'l'); % Thiamine (B1)
model = changeRxnBounds(model, 'EX_ribflv_e', -0.001, 'l'); % Riboflavin (B2)
model = changeRxnBounds(model, 'EX_nac_e', -0.001, 'l'); % Niacin (B3)

% Perform FBA to determine the optimal growth rate
FBAsolution = optimizeCbModel(model, 'max');
if isempty(FBAsolution.f)
    error('The model is infeasible with the given constraints.');
else
    optimalGrowthRate = FBAsolution.f;
    disp(['Optimal growth rate: ', num2str(optimalGrowthRate)]);
end

% Perform FVA to determine the range of fluxes
[minFlux, maxFlux] = fluxVariability(model);

% Display a summary of FVA results
disp('Flux Variability Analysis Results (first 10 reactions):');
for i = 1:min(10, length(model.rxns))
    disp(['Reaction: ', model.rxns{i}, ' Min Flux: ', num2str(minFlux(i)), ' Max Flux: ', num2str(maxFlux(i))]);
end

% Check the feasibility of the model by ensuring no reaction has infeasible fluxes
if any(isnan(minFlux)) || any(isnan(maxFlux))
    error('There are infeasible reactions in the model with the given constraints.');
else
    disp('The model is feasible and ready for sampling.');
end

% Add loop law constraints to the model
model = addLoopLawConstraints(model);

% Perform loopless sampling
numSamples = 1000; % Number of samples
[sampleStruct, mixedFraction] = sampleCbModel(model, [], numSamples, 'loopless', true);

% Display the first few samples
disp('Sampling Results (first 5 samples):');
disp(sampleStruct.samples(:, 1:5));
