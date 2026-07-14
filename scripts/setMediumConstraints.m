clear; close all; clc;

% load('driverModel_WithoutO2','model_7_ctrl', 'model_7_sc1', 'model_7_sc12', 'model_7_sc2',...
%                              'model_H_ctrl', 'model_H_sc1', 'model_H_sc12', 'model_H_sc2');
                         
load('driverModel_WithoutO2onoffLIHC','model_LIHC_ON', 'model_LIHC_OFF');                       
load('driverModel_WithoutO2onoffKIRC','model_KIRC_ON', 'model_KIRC_OFF');  
load('driverModel_WithoutO2onoffKICH','model_KICH_ON', 'model_KICH_OFF');  
load('driverModel_WithoutO2onoffKIRP','model_KIRP_ON', 'model_KIRP_OFF'); 

% Define the list of model names for type 7 models and H models

% models_7 = {model_7_ctrl, model_7_sc1, model_7_sc12, model_7_sc2};
% models_H = {model_H_ctrl, model_H_sc1, model_H_sc12, model_H_sc2};
% 
% modelNames_7 = {'model_7_ctrl', 'model_7_sc1', 'model_7_sc12', 'model_7_sc2'};
% modelNames_H = {'model_H_ctrl', 'model_H_sc1', 'model_H_sc12', 'model_H_sc2'};

models_Kidney = {model_KIRC_ON, model_KIRC_OFF, model_KIRP_ON, model_KIRP_OFF, model_KICH_ON, model_KICH_OFF};
models_Liver = {model_LIHC_ON, model_LIHC_OFF};

modelNames_Kidney = {'model_KIRC_ON', 'model_KIRC_OFF', 'model_KIRP_ON', 'model_KIRP_OFF', 'model_KICH_ON', 'model_KICH_OFF'};
modelNames_Liver  = {'model_LIHC_ON', 'model_LIHC_OFF'};

model_constrained=struct();

%%
medium_composition_RPMI = {'EX_gly[e]';'EX_ala_L[e]';'EX_arg_L[e]';'EX_asn_L[e]';'EX_asp_L[e]';...
                           'EX_cys_L[e]';'EX_glu_L[e]';'EX_his_L[e]';'EX_4HPRO[e]';'EX_ile_L[e]';...
                           'EX_leu_L[e]';'EX_lys_L[e]';'EX_met_L[e]';'EX_phe_L[e]';'EX_pro_L[e]';...
                           'EX_ser_L[e]';'EX_thr_L[e]';'EX_trp_L[e]';'EX_tyr_L[e]';'EX_val_L[e]';...
                           'EX_btn[e]';'EX_chol[e]';'EX_pnto_R[e]';'EX_fol[e]';'EX_ncam[e]';...
                           'EX_pydxn[e]';'EX_ribflv[e]';'EX_thm[e]';'EX_inost[e]';'EX_ca2[e]';...
                           'EX_so4[e]';'EX_k[e]';'EX_hco3[e]';'EX_na1[e]';'EX_pi[e]';...
                           'EX_glc_D[e]';'EX_gln_L[e]';'EX_gthrd[e]'};
                       
met_Conc_mM_RPMI = [0.13333334;2.0552995;1.1494253;0.37878788;0.15037593;...
                        0.20833333;0.13605443;0.09677419;0.15267175;0.3816794;...
                        0.3816794;0.21857923;0.10067114;0.09090909;0.17391305;...
                        0.2857143;0.16806723;0.024509804;0.110497236;0.17094018;...
                        8.20E-04;0.021428572;0.000524;0.002267574;0.008196721;...
                        0.004854369;0.000532;0.002967359;0.19444445;0.42372882;...
                        0.40650406;5.3333335;23.809525;103.44827;5.633803;...
                        11.111111;2.0552995;0.0032573289];
                                        
 medium_composition_DMEM = {'EX_gly[e]';'EX_ala_L[e]';'EX_arg_L[e]';'EX_cys_L[e]';'EX_his_L[e]';...
                                'EX_ile_L[e]';'EX_leu_L[e]';'EX_lys_L[e]';'EX_met_L[e]';'EX_phe_L[e]';...
                                'EX_ser_L[e]';'EX_thr_L[e]';'EX_trp_L[e]';'EX_tyr_L[e]';'EX_val_L[e]';...
                                'EX_chol[e]';'EX_pnto_R[e]';'EX_fol[e]';'EX_ncam[e]';'EX_pydxn[e]';...
                                'EX_ribflv[e]';'EX_thm[e]';'EX_inost[e]';'EX_ca2[e]';'EX_fe3[e]';...
                                'EX_so4[e]';'EX_k[e]';'EX_hco3[e]';'EX_na1[e]';'EX_pi[e]';...
                                'EX_glc_D[e]';'EX_pyr[e]';'EX_gln_L[e]'};
                            
 met_Conc_mM_DMEM = [0.4;3.9723501;0.39810428;0.20127796;0.2
                         0.8015267;0.8015267;0.7978142;0.20134228;0.4
                         0.4;0.79831934;0.078431375;0.39846742;0.8034188
                         0.028571429;0.008385744;0.009070295;0.032786883;0.019417476
                         0.00106383;0.011869436;0.04;1.8018018;2.48E-04
                         0.8139166;5.3333335;44.04762;110.344826;0.9057971;
                         25;1;3.9723501];
                    
%% Loop over the medium compositions for RPMI in type 7 models
for ii = 1:numel(medium_composition_RPMI)
    rxn = medium_composition_RPMI{ii};
    conc = -met_Conc_mM_RPMI(ii);
    
    for i = 1:numel(models_Kidney)
        model = models_Kidney{i};
        modelName = modelNames_Kidney{i};
        %modelSampling = samplingResults(i).modelSampling;
        
        % Find the index of the reaction corresponding to the metabolite in the model
        rxnIndex = find(strcmp(model.rxns, rxn));
        if ~isempty(rxnIndex)
            model.lb(rxnIndex) =  conc;
            disp(model.lb(ismember(model.rxns,'EX_glc_D[e]')));
        else
            disp(['Reaction not found in ' modelName ': ' rxn]);
        end
        
        models_Kidney_RPMI{i}= model;% save here every new model transformed with the addition of the media
        
    end
end

save models_Kidney_RPMI models_Kidney_RPMI

%% Loop over the medium compositions for DMEN in H models
for ii = 1:numel(medium_composition_DMEM)
    rxn = medium_composition_DMEM{ii};
    conc = -met_Conc_mM_DMEM(ii);
    
    for i = 1:numel(models_Liver)
        modelName = modelNames_Liver{i};
             model = models_Liver{i};
        
        % Find the index of the reaction corresponding to the metabolite in the model
        rxnIndex = find(strcmp(model.rxns, rxn));
        if ~isempty(rxnIndex)
           model.lb(rxnIndex) = conc;
            disp(model.lb(ismember(model.rxns,'EX_glc_D[e]'))); 
            disp(model.lb(ismember(model.rxns,'EX_gln_L[e]')));
        else
            disp(['Reaction not found in ' modelName ': ' rxn]);
        end
        
        models_Liver_DMEM{i}=model;
    end
end

save models_Liver_DMEM models_Liver_DMEM

%% Huh H
ResultsLiver=struct();TFluxlLiver=struct();
ResultsKidney=struct();TFluxKidney=struct();

for i=1:numel(models_Liver_DMEM)
changeCobraSolver('glpk');
model=models_Liver_DMEM{i};
model.lb(ismember(model.rxns,'EX_HC02161[e]'))=0;
model.lb(ismember(model.rxns,'EX_peplys[e]'))=0;
%model.lb(ismember(model.rxns,'EX_hdl_hs[e]'))=0;
%model.lb(ismember(model.rxns,'EX_ldl_hs[e]'))=0;
model.lb(ismember(model.rxns,'EX_o2s[e]'))=0;
model.lb(ismember(model.rxns,'EX_o2[e]'))=-50;%Glc x2
model.lb(ismember(model.rxns,'EX_co2[e]'))=-50; %Glc x2

changeObjective(model,'biomass_reaction');
FBALiver=optimizeCbModel(model, 'max', 'zero');
FBALiver.f
[minFluxLiver, maxFluxLiver] = fluxVariability(model, 90);
TFlux=table(model.rxns,FBALiver.x,minFluxLiver,maxFluxLiver);
TFlux.(modelNames_Liver{i})=TFlux;
TFluxEx=TFlux(contains(model.rxns, 'EX_'),:);
TFluxEx = sortrows(TFluxEx,'Var2','ascend');
ResultsLiver.(modelNames_Liver{i})=TFluxEx;
models_Liver_DMEMbiomass{i}=model;

end
save models_Liver_DMEMbiomass models_Liver_DMEMbiomass
%% 7690 7 
for i=1:numel(models_Kidney_RPMI)
changeCobraSolver('glpk')
model=models_Kidney_RPMI{i};
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
models_Kidney_RPMIbiomass{i}=model;
end

save models_Kidney_RPMIbiomass models_Kidney_RPMIbiomass
%%
save setMediumConstraintsONOFF_TCGAsamples.mat
%% 
%Define the list of exchange reactions to include
desiredReactions = {'EX_hco3[e]', 'EX_glc_D[e]', 'EX_o2[e]', 'EX_pyr[e]', 'EX_gln_L[e]', 'EX_co2[e]', 'EX_lac_L[e]', 'EX_h2o[e]', 'EX_ac[e]'};

%% 
modelsFluxEX=struct();
% Create histograms for "7" models
for i = 1:numel(models_Kidney_RPMIbiomass)
    model = models_Kidney_RPMIbiomass{i};
    modelName = modelNames_Kidney{i}; %modelNames7str
    
    % Extract exchange reactions fluxes
    TFluxEx = ResultsKidney.(modelName);
    
    % Filter out non-zero exchange reactions that are in the desiredReactions list
    nonZeroTFluxEx = TFluxEx(ismember(TFluxEx.Var1, desiredReactions) & TFluxEx.Var2 ~= 0, :);
    
    % Check if there are non-zero exchange reactions
    if ~isempty(nonZeroTFluxEx)
        
        modelsFluxEX(i).flux=nonZeroTFluxEx;
        modelsFluxEX(i).modelname=modelName;
        % Create a new figure for the histogram
    end
    
    
    figure;
    
    % Create a histogram of non-zero exchange reaction fluxes
    
    %bar(nonZeroTFluxEx.Var2); % Adjust the number of bins as needed
    
    bar(modelsFluxEX(i).flux.Var2,"yellow")%,modelsFluxEX(2).flux.Var2,"blue",modelsFluxEX(3).flux.Var2,"green",modelsFluxEX(4).flux.Var2,"yelow");
    
    title(['Exchange Reaction Fluxes Histogram for ' modelsFluxEX(i).modelname])%,modelsFluxEX(2).modelname,modelsFluxEX(3).modelname,modelsFluxEX(4).modelname]);
    xlabel('Flux');
    ylabel('Frequency');
    
    % Set X-axis labels to the corresponding reaction IDs
    ax = gca;
    ax.XTickLabel = modelsFluxEX(i).flux.Var1;
    
    % Rotate X-axis labels for better readability
    xtickangle(90);
    
    % Save the histogram as an image
    saveas(gcf, [modelsFluxEX(i).modelname '_EX_histogram.png']);
    saveas(gcf, [modelsFluxEX(i).modelname '_EX_histogram.fig']);
end
%% 
% Create histograms for "H" models
for i = 1:numel(models_Liver_DMEMbiomass)
    model = models_Liver_DMEMbiomass{i};
    modelName = modelNames_Liver{i};

    % Extract exchange reactions fluxes
    TFluxEx = ResultsLiver.(modelName);

    % Filter out non-zero exchange reactions that are in the desiredReactions list
    nonZeroTFluxEx = TFluxEx(ismember(TFluxEx.Var1, desiredReactions) & TFluxEx.Var2 ~= 0, :);

    % Check if there are non-zero exchange reactions
    if ~isempty(nonZeroTFluxEx)
        % Create a new figure for the histogram
        figure;

        % Create a histogram of non-zero exchange reaction fluxes
        bar(nonZeroTFluxEx.Var2); % Adjust the number of bins as needed
        title(['Non-Zero Exchange Reaction Fluxes Histogram for ' modelName]);
        xlabel('Flux');
        ylabel('Frequency');

        % Set X-axis labels to the corresponding reaction IDs
        ax = gca;
        ax.XTickLabel = nonZeroTFluxEx.Var1;

        % Rotate X-axis labels for better readability
        xtickangle(90);

         % Save the histogram as an image
        saveas(gcf, [modelName '_nonzero_histogram.png']);
        saveas(gcf, [modelName '_nonzero_histogram.fig']);
    end
end

