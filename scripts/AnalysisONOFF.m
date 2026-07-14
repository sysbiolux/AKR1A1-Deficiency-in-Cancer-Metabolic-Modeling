clear; close all; clc;
feature astheightlimit 2000
solverOK=1;

% Main script for metabolic modeling analysis

%%
% Load or generate metabolic network model
load ('driverModel_WithoutO2','colnamesH','colnames7','consistent_model','models_keep_7sample','models_keep_Hsample');

%load modelNames_Hwithouto2s.mat
%load modelNames_7withouto2s.mat
load GeneDrugRelationsUpdate.mat
load models_H_DMEMbiomass.mat
load models_7_RPMIbiomass.mat

model_7_ctrl = models_7_RPMIbiomass{1,1};
model_7_sc1  = models_7_RPMIbiomass{1,2}; 
model_7_sc12 = models_7_RPMIbiomass{1,3};
model_7_sc2  = models_7_RPMIbiomass{1,4};   
model_H_ctrl = models_H_DMEMbiomass{1,1};
model_H_sc1 = models_H_DMEMbiomass{1,2};
model_H_sc12 = models_H_DMEMbiomass{1,3};
model_H_sc2 = models_H_DMEMbiomass{1,4}; 

modelNames_7 = {'model_7_ctrl', 'model_7_sc1', 'model_7_sc12', 'model_7_sc2'}; 
modelNames_H = {'model_H_ctrl', 'model_H_sc1', 'model_H_sc12', 'model_H_sc2'};

% Load or generate metabolic network model
%load ('driverModel_WithoutO2','colnamesH','colnames7','consistent_model','models_keep_7sample','models_keep_Hsample');
%load ('driverModel_WithoutO2offALL');
%load modelNames_Hwithouto2s.mat
%load modelNames_7withouto2s.mat

%load GeneDrugRelations.mat
load GeneDrugRelationsUpdate.mat

% load models_H_DMEMbiomass.mat
% load models_7_RPMIbiomass.mat

load models_Kidney_RPMIbiomass
load models_Liver_DMEMbiomass

model_KIRC_ON = models_Kidney_RPMIbiomass{1,1};
model_KIRC_OFF  = models_Kidney_RPMIbiomass{1,2}; 
model_KIRP_ON = models_Kidney_RPMIbiomass{1,3};
model_KIRP_OFF  = models_Kidney_RPMIbiomass{1,4};   
model_KICH_ON = models_Kidney_RPMIbiomass{1,5};
model_KICH_OFF = models_Kidney_RPMIbiomass{1,6};
model_LIHC_ON = models_Liver_DMEMbiomass{1,1};
model_LIHC_OFF = models_Liver_DMEMbiomass{1,2};

modelNames_kidney = {'model_KIRC_ON', 'model_KIRC_OFF', 'model_KIRP_ON', 'model_KIRP_OFF', 'model_KICH_ON', 'model_KICH_OFF'}; 
modelNames_liver = {'model_LIHC_ON', 'model_LIHC_OFF'};

load discretizedOFF_KICH
load discretizedON_KICH
load discretizedOFF_KIRC
load discretizedON_KIRC
load discretizedOFF_KIRP
load discretizedON_KIRP
load discretizedOFF_LIHC
load discretizedON_LIHC

feature astheightlimit 2000

% %% Jaccard Similarity Analysis
% feature astheightlimit 2000
% altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
%     255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar
% 
% %%
% J = squareform(pdist(models_keep_Hsample','jaccard'));
% %Jaccard similarity plots for sample models H
% 
% cgo_J = clustergram(1-J,...
%     'RowLabels', colnamesH,...
%     'ColumnLabels', colnamesH,...
%     'ColumnLabelsRotate',45, ...
%     'Cluster', 'all', ...
%     'symmetric','False',...
%     'Colormap', altcolor);
% addTitle(cgo_J,{'Model similarity based on Jaccard distance','models_keep'})
% 
% % Wait for the clustergram to be created
% drawnow; % Ensure that the clustergram is rendered
% 
% % Get handle to the clustergram figure
% cgf = plot(cgo_J); % This should be a figure handle
% 
% 
% % Save the figure with color scale
% saveas(cgf, 'ModelsimilaritybasedonJaccarddistance_H.png');
% saveas(gcf, 'ModelsimilaritybasedonJaccarddistance_H.pdf');
% 
% %%
% J = squareform(pdist(models_keep_7sample','jaccard'));
% %Jaccard similarity plots for sample models 7
% 
% cgo_J = clustergram(1-J,...
%     'RowLabels', colnames7,...
%     'ColumnLabels', colnames7,...
%     'ColumnLabelsRotate',45, ...
%     'Cluster', 'all', ...
%     'symmetric','False',...
%     'Colormap', altcolor);
% addTitle(cgo_J,{'Model similarity based on Jaccard distance','models_keep'})
% 
% % Wait for the clustergram to be created
% drawnow; % Ensure that the clustergram is rendered
% 
% % Get handle to the clustergram figure
% cgf = plot(cgo_J); % This should be a figure handle
% 
% % Add colorbar to the figure
% colorbar(cgf,'eastoutside');
% 
% % Save the figure with color scale
% saveas(cgf, 'ModelsimilaritybasedonJaccarddistance_7.png');
% saveas(gcf, 'ModelsimilaritybasedonJaccarddistance_h.pdf');
% 
% %% Pathway analysis 
% Pathways = table(unique(consistent_model.subSystems));
% [pathways, ~, ub] = unique(consistent_model.subSystems);
% path_counts = histc(ub, 1:length(pathways));
% T = table(pathways, path_counts);
% [~, ia, ib] = intersect(Pathways.Var1, T.pathways);
% Pathways.consistent(ia) = T.path_counts(ib);
% Pathways.Properties.VariableNames{1}='Pathways';
% 
% %% 
% PathwaysH=Pathways;
% % pathway information for the consensus models for H 
% for i=1: numel(colnamesH)
%    
% [pathways, ~, ub] = unique(consistent_model.subSystems(models_keep_Hsample(:,i)~=0));
% path_counts = histc(ub, 1:length(pathways));
% T = table(pathways, path_counts);
% [~, ia, ib] = intersect(PathwaysH.Pathways, T.pathways);
% PathwaysH.Var2(ia) = T.path_counts(ib) ;
% PathwaysH.Properties.VariableNames{2+i} = colnamesH{i};
% end
% 
% Pathways7=Pathways;
% % pathway information for the consensus models for 7 
% for i=1: numel(colnames7)
% [pathways, ~, ub] = unique(consistent_model.subSystems(models_keep_7sample(:,i)~=0));
% path_counts = histc(ub, 1:length(pathways));
% T = table(pathways, path_counts);
% [~, ia, ib] = intersect(Pathways7.Pathways, T.pathways);
% Pathways7.Var2(ia) = T.path_counts(ib) ;
% Pathways7.Properties.VariableNames{2+i} = colnamesH{i};
% end
% 
% %% plotting pathways
% PathwayActivityH = PathwaysH;
% for i=3:size(PathwayActivityH,2)
%     PathwayActivityH(:,i) = array2table(table2array(PathwayActivityH(:,i))./table2array(PathwayActivityH(:,2)));
% end
% 
% cgo = clustergram(table2array(PathwayActivityH(:,3:end)),...
%     'RowLabels', PathwayActivityH.Pathways,...
%     'ColumnLabels', PathwayActivityH.Properties.VariableNames(3:end),...
%     'ColumnLabelsRotate',270, ...
%     'Cluster', 'all', ...
%     'symmetric','False',...
%     'Colormap', altcolor);
% h = plot(cgo); set(h,'TickLabelInterpreter','none');
% colorbar(h)
% title(h,'Pathway activity for all models H')
% 
% % Wait for the clustergram to be created
% drawnow; % Ensure that the clustergram is rendered
% 
% % Get handle to the clustergram figure
% cgf = plot(cgo); % This should be a figure handle
% 
% % Add colorbar to the figure
% colorbar(cgf,'eastoutside');
% 
% % Save the figure with color scale
% saveas(cgf, 'Pathwayactivityforallmodels_H.png');
% 
% 
% PathwayActivity7 = Pathways7;
% for i=3:size(PathwayActivity7,2)
%     PathwayActivity7(:,i) = array2table(table2array(PathwayActivity7(:,i))./table2array(PathwayActivity7(:,2)));
% end
% 
% cgo = clustergram(table2array(PathwayActivity7(:,3:end)),...
%     'RowLabels', PathwayActivity7.Pathways,...
%     'ColumnLabels', PathwayActivity7.Properties.VariableNames(3:end),...
%     'ColumnLabelsRotate',270, ...
%     'Cluster', 'all', ...
%     'symmetric','False',...
%     'Colormap', altcolor);
% h = plot(cgo); set(h,'TickLabelInterpreter','none');
% colorbar(h)
% title(h,'Pathway activity for all models 7')
% 
% % Wait for the clustergram to be created
% drawnow; % Ensure that the clustergram is rendered
% 
% % Get handle to the clustergram figure
% cgf = plot(cgo); % This should be a figure handle
% 
% % Add colorbar to the figure
% colorbar(cgf,'eastoutside');
% 
% % Save the figure with color scale
% saveas(cgf, 'Pathwayactivityforallmodels_7.png');

%% Single Gene Deletion Analysis
geneList_keep=struct();
%% run the single gene deletion for 769-P

for i=1:numel(modelNames_kidney)
    name =  modelNames_kidney(i);
    model = models_Kidney_RPMIbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_biomass),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRateWT_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('essential', modelNames_kidney(i));
    save (name2{1},'TEssential');

end
save singlegenedeletionKidney
%% run the single gene deletion for Huh7
for i=1:numel(modelNames_liver)
    name =  modelNames_liver(i);
    model = models_Liver_DMEMbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_ATP),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('essential', modelNames_liver(i));
    save (name2{1},'TEssential');

end
save singlegenedeletionLiver

%% Drug deletion
% define a list of drugs
load GeneDrugRelationsUpdate.mat
DrugList = unique(GeneDrugRelations.DrugName);
%% Huh7

model_H_ctrl= changeObjective(model_H_ctrl,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_H_ctrl,'FBA',DrugList);
    
Drug_grRatio_biomass(:,1)    = grRatio;
Drug_grRateKO_biomass(:,1)   = grRateKO;
Drug_grRateWT_biomass(:,1)   = grRateWT;
    
control_H_drug= DrugList(Drug_grRatio_biomass(:,1) ==0);
save control_H_drug control_H_drug
%%
model_H_sc1= changeObjective(model_H_sc1,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_H_sc1,'FBA',DrugList);
    
Drug_grRatio_biomass(:,2)    = grRatio;
Drug_grRateKO_biomass(:,2)   = grRateKO;
Drug_grRateWT_biomass(:,2)   = grRateWT;
    
Sc1_H_Drug = DrugList(Drug_grRatio_biomass(:,2) ==0);
save Sc1_H_Drug Sc1_H_Drug
%%
model_H_sc2= changeObjective(model_H_sc2,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_H_sc2,'FBA',DrugList);
    
Drug_grRatio_biomass(:,3)    = grRatio;
Drug_grRateKO_biomass(:,3)   = grRateKO;
Drug_grRateWT_biomass(:,3)   = grRateWT;

Sc2_H_Drug = DrugList(Drug_grRatio_biomass(:,3) ==0);
save Sc2_H_Drug Sc2_H_Drug
%%
model_H_sc12= changeObjective(model_H_sc12,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_H_sc12,'FBA',DrugList);
    
Drug_grRatio_biomass(:,4)    = grRatio;
Drug_grRateKO_biomass(:,4)   = grRateKO;
Drug_grRateWT_biomass(:,4)   = grRateWT;
    
Sc12_H_drug = DrugList(Drug_grRatio_biomass(:,4) ==0);
save Sc12_H_drug Sc12_H_drug


%% Drug Deletion for 769-P
%load GeneDrugRelations.mat
DrugList = unique(GeneDrugRelations.DrugName);
model_7_ctrl= changeObjective(model_7_ctrl,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_ctrl,'FBA',DrugList);
    
Drug_grRatio_biomass(:,1)    = grRatio;
Drug_grRateKO_biomass(:,1)   = grRateKO;
Drug_grRateWT_biomass(:,1)   = grRateWT;
    
control_7_Drug= DrugList(Drug_grRatio_biomass(:,1) ==0);
save control_7_Drug control_7_Drug

%%
model_7_sc1= changeObjective(model_7_sc1,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc1,'FBA',DrugList);
    
Drug_grRatio_biomass(:,2)    = grRatio;
Drug_grRateKO_biomass(:,2)   = grRateKO;
Drug_grRateWT_biomass(:,2)   = grRateWT;
    
Sc1_7_Drug = DrugList(Drug_grRatio_biomass(:,2) ==0);
save Sc1_7_Drug Sc1_7_Drug
%%
model_7_sc2= changeObjective(model_7_sc2,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc2,'FBA',DrugList);
    
Drug_grRatio_biomass(:,3)    = grRatio;
Drug_grRateKO_biomass(:,3)   = grRateKO;
Drug_grRateWT_biomass(:,3)   = grRateWT;

    
Sc2_7_Drug = DrugList(Drug_grRatio_biomass(:,3) ==0);
save Sc2_7_Drug Sc2_7_Drug

%%
model_7_sc12= changeObjective(model_7_sc12,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc12,'FBA',DrugList);
    
Drug_grRatio_biomass(:,4)    = grRatio;
Drug_grRateKO_biomass(:,4)   = grRateKO;
Drug_grRateWT_biomass(:,4)   = grRateWT;
    
Sc12_7_Drug = DrugList(Drug_grRatio_biomass(:,4) ==0);
save Sc12_7_Drug Sc12_7_Drug 

save workspace_singleDrugDeletionCELLLINES

%% Drug deletion
% define a list of drugs
load GeneDrugRelationsUpdate.mat
DrugList = unique(GeneDrugRelations.DrugName);
%%
model_LIHC_ON= changeObjective(model_LIHC_ON,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_LIHC_ON,'FBA',DrugList);
    
Drug_grRatio_biomass(:,1)    = grRatio;
Drug_grRateKO_biomass(:,1)   = grRateKO;
Drug_grRateWT_biomass(:,1)   = grRateWT;

LIHC_ON_Drug = DrugList(Drug_grRatio_biomass(:,1) ==0);
save LIHC_ON_Drug LIHC_ON_Drug
%%
model_LIHC_OFF= changeObjective(model_LIHC_OFF,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_LIHC_OFF,'FBA',DrugList);
    
Drug_grRatio_biomass(:,2)    = grRatio;
Drug_grRateKO_biomass(:,2)   = grRateKO;
Drug_grRateWT_biomass(:,2)   = grRateWT;
    
LIHC_OFF_Drug = DrugList(Drug_grRatio_biomass(:,2) ==0);
save LIHC_OFF_Drug LIHC_OFF_Drug
%% Drug Deletion for 769-P
%load GeneDrugRelations.mat
DrugList = unique(GeneDrugRelations.DrugName);
%%
model_KIRC_ON= changeObjective(model_KIRC_ON,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KIRC_ON,'FBA',DrugList);
    
Drug_grRatio_biomass(:,1)    = grRatio;
Drug_grRateKO_biomass(:,1)   = grRateKO;
Drug_grRateWT_biomass(:,1)   = grRateWT;
    
KIRC_ON_Drug= DrugList(Drug_grRatio_biomass(:,1) ==0);
save KIRC_ON_Drug KIRC_ON_Drug

%%
model_KIRC_OFF= changeObjective(model_KIRC_OFF,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KIRC_OFF,'FBA',DrugList);
    
Drug_grRatio_biomass(:,2)    = grRatio;
Drug_grRateKO_biomass(:,2)   = grRateKO;
Drug_grRateWT_biomass(:,2)   = grRateWT;
    
    
KIRC_OFF_Drug = DrugList(Drug_grRatio_biomass(:,2) ==0);
save KIRC_OFF_Drug KIRC_OFF_Drug
%%
model_KIRP_ON= changeObjective(model_KIRP_ON,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KIRP_ON,'FBA',DrugList);
    
Drug_grRatio_biomass(:,3)    = grRatio;
Drug_grRateKO_biomass(:,3)   = grRateKO;
Drug_grRateWT_biomass(:,3)   = grRateWT;

    
KIRP_ON_Drug = DrugList(Drug_grRatio_biomass(:,3) ==0);
save KIRP_ON_Drug KIRP_ON_Drug
%%
model_KIRP_OFF= changeObjective(model_KIRP_OFF,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KIRP_OFF,'FBA',DrugList);
    
Drug_grRatio_biomass(:,4)    = grRatio;
Drug_grRateKO_biomass(:,4)   = grRateKO;
Drug_grRateWT_biomass(:,4)   = grRateWT;

    
KIRP_OFF_Drug = DrugList(Drug_grRatio_biomass(:,4) ==0);
save KIRP_OFF_Drug KIRP_OFF_Drug
%%
model_KICH_ON= changeObjective(model_KICH_ON,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KICH_ON,'FBA',DrugList);
    
Drug_grRatio_biomass(:,5)    = grRatio;
Drug_grRateKO_biomass(:,5)   = grRateKO;
Drug_grRateWT_biomass(:,5)   = grRateWT;

    
KICH_ON_Drug = DrugList(Drug_grRatio_biomass(:,5) ==0);
save KICH_ON_Drug KICH_ON_Drug
%%
model_KICH_OFF= changeObjective(model_KICH_OFF,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_KICH_OFF,'FBA',DrugList);
    
Drug_grRatio_biomass(:,6)    = grRatio;
Drug_grRateKO_biomass(:,6)   = grRateKO;
Drug_grRateWT_biomass(:,6)   = grRateWT;

    
KICH_OFF_Drug = DrugList(Drug_grRatio_biomass(:,6) ==0);
save KICH_OFF_Drug KICH_OFF_Drug

save workspace_singleDrugDeletionTUMORSAMPLES
%% Similarity Based on Flux Variability Analysis
%FVA 7
[minFlux_kichon, maxFlux_kichon] = fluxVariability(model_KICH_ON);
[minFlux_kichoff, maxFlux_kichoff] = fluxVariability(model_KICH_OFF);   
[minFlux_kircon, maxFlux_kircon] = fluxVariability(model_KIRC_ON); 
[minFlux_kircoff, maxFlux_kircoff] = fluxVariability(model_KIRC_OFF);    
[minFlux_kirpon, maxFlux_kirpon] = fluxVariability(model_KIRP_ON); 
[minFlux_kirpoff, maxFlux_kirpoff] = fluxVariability(model_KIRP_OFF); 

T_fluxes_kichon=table(model_KICH_ON.rxns,model_KICH_ON.rxnNames,minFlux_kichon,maxFlux_kichon);
T_fluxes_kichoff=table(model_KICH_OFF.rxns,model_KICH_OFF.rxnNames,minFlux_kichoff,maxFlux_kichoff);
T_fluxes_kircon=table(model_KIRC_ON.rxns,model_KIRC_ON.rxnNames,maxFlux_kircon);
T_fluxes_kircoff=table(model_KIRC_OFF.rxns,model_KIRC_OFF.rxnNames,maxFlux_kircoff);
T_fluxes_kirpon=table(model_KIRP_ON.rxns,model_KIRP_ON.rxnNames,maxFlux_kirpon);
T_fluxes_kirpoff=table(model_KIRP_OFF.rxns,model_KIRP_OFF.rxnNames,maxFlux_kirpoff);

[minFlux_lihcon, maxFlux_lihcon] = fluxVariability(model_LIHC_ON);
[minFlux_lihcoff, maxFlux_lihcoff] = fluxVariability(model_LIHC_OFF);

T_fluxes_lihcon=table(model_LIHC_ON.rxns,model_LIHC_ON.rxnNames,minFlux_lihcon,maxFlux_lihcon);
T_fluxes_lihcoff=table(model_LIHC_OFF.rxns,model_LIHC_OFF.rxnNames,minFlux_lihcoff,maxFlux_lihcoff);


save analysis_part1onoff
%%
% Initialize matrices to store min and max flux values
minFlux_keep = zeros(numel(consistent_model.rxns), 6);
maxFlux_keep = zeros(numel(consistent_model.rxns), 6);

% Store the FVA results in the matrices for ON and OFF models
minFlux_keep(ismember(consistent_model.rxns, model_KICH_ON.rxns), 1) = minFlux_kichon;
maxFlux_keep(ismember(consistent_model.rxns, model_KICH_ON.rxns), 1) = maxFlux_kichon;
minFlux_keep(ismember(consistent_model.rxns, model_KICH_OFF.rxns), 2) = minFlux_kichoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KICH_OFF.rxns), 2) = maxFlux_kichoff;

minFlux_keep(ismember(consistent_model.rxns, model_KIRC_ON.rxns), 3) = minFlux_kircon;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRC_ON.rxns), 3) = maxFlux_kircon;
minFlux_keep(ismember(consistent_model.rxns, model_KIRC_OFF.rxns), 4) = minFlux_kircoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRC_OFF.rxns), 4) = maxFlux_kircoff;

minFlux_keep(ismember(consistent_model.rxns, model_KIRP_ON.rxns), 5) = minFlux_kirpon;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRP_ON.rxns), 5) = maxFlux_kirpon;
minFlux_keep(ismember(consistent_model.rxns, model_KIRP_OFF.rxns), 6) = minFlux_kirpoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRP_OFF.rxns), 6) = maxFlux_kirpoff;
%%
minFlux_keep(ismember(consistent_model.rxns, model_LIHC_ON.rxns), 7) = minFlux_lihcon;
maxFlux_keep(ismember(consistent_model.rxns, model_LIHC_ON.rxns), 7) = maxFlux_lihcon;
minFlux_keep(ismember(consistent_model.rxns, model_LIHC_OFF.rxns), 8) = minFlux_lihcoff;
maxFlux_keep(ismember(consistent_model.rxns, model_LIHC_OFF.rxns), 8) = maxFlux_lihcoff;
%%
% Initialize matrix to store the similarity results
res_keep = zeros(size(minFlux_keep, 2), size(minFlux_keep, 2));

% Calculate similarity for each pair of models
for i = 1:size(minFlux_keep, 2)
    for j = 1:size(minFlux_keep, 2)
        v1mins = minFlux_keep(:, i);
        v2mins = minFlux_keep(:, j);
        v1maxs = maxFlux_keep(:, i);
        v2maxs = maxFlux_keep(:, j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res_keep(i, j) = res;
    end
end

% Initialize matrix to store subsystem-specific similarity results
res2_keep = zeros(numel(uniSys), size(minFlux_keep, 2) * (size(minFlux_keep, 2) - 1) / 2);
k = 0;

% Calculate subsystem-specific similarity for each pair of models
for i = 1:size(minFlux_keep, 2)
    for j = 1:size(minFlux_keep, 2)
        if i > j
            k = k + 1;
            for counter = 1:numel(uniSys)
                match = ismember(consistent_model.subSystems, uniSys(counter));
                v1mins = minFlux_keep(match, i);
                v2mins = minFlux_keep(match, j);
                v1maxs = maxFlux_keep(match, i);
                v2maxs = maxFlux_keep(match, j);
                res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
                res2_keep(counter, k) = res;
            end
        end
    end
end

% Create and save the table with subsystem similarity results
T_FVA_ONOFF = table(uniSys, res2_keep);
save('T_FVA_ONOFF.mat', 'T_FVA_ONOFF');

%%
minFlux_7_keep=zeros(numel(consistent_model.rxns),6);
maxFlux_7_keep=zeros(numel(consistent_model.rxns),6);
minFlux_7_keep(ismember(consistent_model.rxns,model_KICH_ON.rxns),1)=minFlux_kichon;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KICH_ON.rxns),1)=maxFlux_kichon;
minFlux_7_keep(ismember(consistent_model.rxns,model_KICH_OFF.rxns),2)=minFlux_kichoff;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KICH_OFF.rxns),2)=maxFlux_kichoff;
minFlux_7_keep(ismember(consistent_model.rxns,model_KIRC_ON.rxns),3)=minFlux_kircon;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KIRC_ON.rxns),3)=maxFlux_kircon;
minFlux_7_keep(ismember(consistent_model.rxns,model_KIRC_OFF.rxns),4)=minFlux_kircoff;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KIRC_OFF.rxns),4)=maxFlux_kircoff;
minFlux_7_keep(ismember(consistent_model.rxns,model_KIRP_ON.rxns),5)=minFlux_kirpon;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KIRP_ON.rxns),5)=maxFlux_kirpon;
minFlux_7_keep(ismember(consistent_model.rxns,model_KIRP_OFF.rxns),6)=minFlux_kirpoff;
maxFlux_7_keep(ismember(consistent_model.rxns,model_KIRP_OFF.rxns),6)=maxFlux_kirpoff;
res_keep=zeros(size(minFlux_7_keep,2),size(minFlux_7_keep,2));
for i=1:size(minFlux_7_keep,2)
    for j=1:size(minFlux_7_keep,2)
        v1mins=minFlux_7_keep(:,i);
        v2mins=minFlux_7_keep(:,j);
        v1maxs=maxFlux_7_keep(:,i);
        v2maxs=maxFlux_7_keep(:,j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res_keep(i,j)=res;
    end
end
k=0;
res2_keep=zeros(1,3);
%disp('contol-s1', 'contol-s2', 's1-s2')
uniSys=unique(consistent_model.subSystems);
for i=1:size(minFlux_7_keep,2)
    for j=1:size(minFlux_7_keep,2)
        if i>j
            k=k+1;
            for counter=1:numel(uniSys)
                match=ismember(consistent_model.subSystems,uniSys(counter));
        v1mins=minFlux_7_keep(match,i);
        v2mins=minFlux_7_keep(match,j);
        v1maxs=maxFlux_7_keep(match,i);
        v2maxs=maxFlux_7_keep(match,j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res2_keep(counter,k)=res;
            end
        end
    end
end
T_FVA_7_ONOFF=table(uniSys,res2_keep);

save('T_FVA_7_ONOFF.mat','T_FVA_7_ONOFF')  % function form
%%
%FVA huh7
minFlux_H_keep=zeros(numel(consistent_model.rxns),2);
maxFlux_H_keep=zeros(numel(consistent_model.rxns),2);
minFlux_H_keep(ismember(consistent_model.rxns,model_LIHC_ON.rxns),1)=minFlux_lihcon;
maxFlux_H_keep(ismember(consistent_model.rxns,model_LIHC_ON.rxns),1)=maxFlux_lihcon;
minFlux_H_keep(ismember(consistent_model.rxns,model_LIHC_OFF.rxns),2)=minFlux_lihcoff;
maxFlux_H_keep(ismember(consistent_model.rxns,model_LIHC_OFF.rxns),2)=maxFlux_lihcoff;

res_keep=zeros(size(minFlux_H_keep,2),size(minFlux_H_keep,2));
for i=1:size(minFlux_H_keep,2)
    for j=1:size(minFlux_H_keep,2)
        v1mins=minFlux_H_keep(:,i);
        v2mins=minFlux_H_keep(:,j);
        v1maxs=maxFlux_H_keep(:,i);
        v2maxs=maxFlux_H_keep(:,j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res_keep(i,j)=res;
    end
end
k=0;
res2_keep=zeros(1,2);
%disp('contol-s1'; 'contol-s2'; 's1-s2');
uniSys=unique(consistent_model.subSystems);
for i=1:size(minFlux_H_keep,2)
    for j=1:size(minFlux_H_keep,2)
        if i>j
            k=k+1;
            for counter=1:numel(uniSys)
                match=ismember(consistent_model.subSystems,uniSys(counter));
        v1mins=minFlux_H_keep(match,i);
        v2mins=minFlux_H_keep(match,j);
        v1maxs=maxFlux_H_keep(match,i);
        v2maxs=maxFlux_H_keep(match,j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res2_keep(counter,k)=res;
            end
        end
    end
end
T_FVA_H_ONOFF=table(uniSys,res2_keep);

save('T_FVA_H_ONOFF.mat','T_FVA_H_ONOFF')  

%% heatmap_FVAsimilarity
% Load the similarity results
load('T_FVA_ONOFF.mat', 'T_FVA_ONOFF');

% Check the number of available subsystems
numSubSystems = numel(T_FVA_ONOFF.uniSys);

% Determine the number of pathways to plot
numPathways = min(41, numSubSystems);

% Ensure that res2_keep is available and valid
if numel(T_FVA_ONOFF.res2_keep) < numPathways
    error('Not enough elements in res2_keep to extract the top %d pathways.', numPathways);
end

% Sort the similarity values and get the top indices
[sortedValues, sortIndex] = sort(T_FVA_ONOFF.res2_keep(:), 'descend');

% Filter valid indices
validIndices = sortIndex(sortIndex <= numSubSystems);
topIndex = validIndices(1:numPathways);

% Debugging information to verify sizes and indices
disp(['Number of Subsystems: ', num2str(numSubSystems)]);
disp(['Number of Pathways: ', num2str(numPathways)]);
disp(['Top Index: ', num2str(topIndex')]);

% Get the corresponding subsystems and values
topSubsystems = T_FVA_ONOFF.uniSys(topIndex);
topValues = T_FVA_ONOFF.res2_keep(topIndex);

% Confirm the number of comparisons (pairs of ON/OFF states)
numComparisons = size(T_FVA_ONOFF, 2) - 1; % Adjust based on actual data structure

% Ensure that the length of topValues is compatible with numComparisons
if mod(numel(topValues), numComparisons) ~= 0
    error('The length of topValues is not compatible with the number of comparisons.');
end

% Reshape topValues into a matrix with numComparisons columns
heatmapData = reshape(topValues, numPathways, numComparisons);

% Customize column labels
modelNames = {'KICH_ON_OFF', 'KIRC_ON_OFF', 'KIRP_ON_OFF', 'LIHC_ON_OFF', 'Other1', 'Other2'}; % Adjust as needed

% Verify the dimensions
disp(['Size of heatmapData: ', num2str(size(heatmapData))]);
disp(['Number of model names: ', num2str(numel(modelNames))]);

% Custom colormap
custom_colormap = redbluecmap(101);

% Generate the heatmap
figure;
h = heatmap(modelNames(1:numComparisons), topSubsystems, heatmapData, ...
    'Colormap', custom_colormap, ...
    'ColorLimits', [min(topValues), max(topValues)], ...
    'ColorbarVisible', 'on');

title('Top 41 Different Pathways based on FVA Similarity');
xlabel('Cancer Models');
ylabel('Pathways');

% Save the heatmap
saveas(gcf, 'FVA_similarity_heatmap_top41.png');
saveas(gcf, 'FVA_similarity_heatmap_top41.pdf');
%%
% Load consistent model and previous results
load('consistent_model.mat'); % Assuming the model is stored in a .mat file
load('T_FVA_ONOFF.mat', 'T_FVA_ONOFF'); % Load previous results if necessary

% Initialize matrices to store min and max flux values
minFlux_keep = zeros(numel(consistent_model.rxns), 8);
maxFlux_keep = zeros(numel(consistent_model.rxns), 8);

% Store the FVA results in the matrices for ON and OFF models
minFlux_keep(ismember(consistent_model.rxns, model_KICH_ON.rxns), 1) = minFlux_kichon;
maxFlux_keep(ismember(consistent_model.rxns, model_KICH_ON.rxns), 1) = maxFlux_kichon;
minFlux_keep(ismember(consistent_model.rxns, model_KICH_OFF.rxns), 2) = minFlux_kichoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KICH_OFF.rxns), 2) = maxFlux_kichoff;

minFlux_keep(ismember(consistent_model.rxns, model_KIRC_ON.rxns), 3) = minFlux_kircon;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRC_ON.rxns), 3) = maxFlux_kircon;
minFlux_keep(ismember(consistent_model.rxns, model_KIRC_OFF.rxns), 4) = minFlux_kircoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRC_OFF.rxns), 4) = maxFlux_kircoff;

minFlux_keep(ismember(consistent_model.rxns, model_KIRP_ON.rxns), 5) = minFlux_kirpon;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRP_ON.rxns), 5) = maxFlux_kirpon;
minFlux_keep(ismember(consistent_model.rxns, model_KIRP_OFF.rxns), 6) = minFlux_kirpoff;
maxFlux_keep(ismember(consistent_model.rxns, model_KIRP_OFF.rxns), 6) = maxFlux_kirpoff;

minFlux_keep(ismember(consistent_model.rxns, model_LIHC_ON.rxns), 7) = minFlux_lihcon;
maxFlux_keep(ismember(consistent_model.rxns, model_LIHC_ON.rxns), 7) = maxFlux_lihcon;
minFlux_keep(ismember(consistent_model.rxns, model_LIHC_OFF.rxns), 8) = minFlux_lihcoff;
maxFlux_keep(ismember(consistent_model.rxns, model_LIHC_OFF.rxns), 8) = maxFlux_lihcoff;

% Initialize matrix to store the similarity results
res_keep = zeros(size(minFlux_keep, 2), size(minFlux_keep, 2));

% Calculate similarity for each pair of models
for i = 1:size(minFlux_keep, 2)
    for j = 1:size(minFlux_keep, 2)
        v1mins = minFlux_keep(:, i);
        v2mins = minFlux_keep(:, j);
        v1maxs = maxFlux_keep(:, i);
        v2maxs = maxFlux_keep(:, j);
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res_keep(i, j) = res;
    end
end

% Initialize matrix to store subsystem-specific similarity results
res2_keep = zeros(numel(uniSys), size(minFlux_keep, 2) * (size(minFlux_keep, 2) - 1) / 2);
k = 0;

% Calculate subsystem-specific similarity for each pair of models
for i = 1:size(minFlux_keep, 2)
    for j = 1:size(minFlux_keep, 2)
        if i > j
            k = k + 1;
            for counter = 1:numel(uniSys)
                match = ismember(consistent_model.subSystems, uniSys(counter));
                v1mins = minFlux_keep(match, i);
                v2mins = minFlux_keep(match, j);
                v1maxs = maxFlux_keep(match, i);
                v2maxs = maxFlux_keep(match, j);
                res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
                res2_keep(counter, k) = res;
            end
        end
    end
end

% Create and save the table with subsystem similarity results
T_FVA_ONOFF = table(uniSys, res2_keep);
save('T_FVA_ONOFF.mat', 'T_FVA_ONOFF');

%%
% Load the similarity results
load('T_FVA_ONOFF.mat', 'T_FVA_ONOFF');

% Determine the number of pathways to plot
numPathways = 41;

% Initialize cell array to store the top indices for each comparison
topIndices = cell(4, 1);

% Comparisons: KICH (1,2), KIRC (3,4), KIRP (5,6), LIHC (7,8)
comparisons = [1, 2; 3, 4; 5, 6; 7, 8];

% Identify the top 41 most different pathways for each comparison
for idx = 1:size(comparisons, 1)
    comp1 = comparisons(idx, 1);
    comp2 = comparisons(idx, 2);
    diffs = abs(res2_keep(:, comp1) - res2_keep(:, comp2));
    [~, sortIndex] = sort(diffs, 'descend');
    topIndices{idx} = sortIndex(1:numPathways);
end

% Debugging information to verify sizes and indices
disp('Top Indices for each comparison:');
disp(topIndices);
%%
% Load the similarity results
load('T_FVA_ONOFF.mat', 'T_FVA_ONOFF');

% Determine the number of pathways to plot
numPathways = 41;

% Initialize cell array to store the top indices for each comparison
topIndices = cell(4, 1);

% Comparisons: KICH (1,2), KIRC (3,4), KIRP (5,6), LIHC (7,8)
comparisons = [1, 2; 3, 4; 5, 6; 7, 8];

% Identify the top 41 most different pathways for each comparison
for idx = 1:size(comparisons, 1)
    comp1 = comparisons(idx, 1);
    comp2 = comparisons(idx, 2);
    diffs = abs(res2_keep(:, comp1) - res2_keep(:, comp2));
    [~, sortIndex] = sort(diffs, 'descend');
    topIndices{idx} = sortIndex(1:numPathways);
end

% Debugging information to verify sizes and indices
disp('Top Indices for each comparison:');
disp(topIndices);
%%
% Custom colormap
custom_colormap = redbluecmap(101);

% Model names
modelNames = {'KICH_ON_OFF', 'KIRC_ON_OFF', 'KIRP_ON_OFF', 'LIHC_ON_OFF'};

% Generate heatmaps for each comparison
for idx = 1:size(comparisons, 1)
    comp1 = comparisons(idx, 1);
    comp2 = comparisons(idx, 2);
    topIndex = topIndices{idx};

    % Get the corresponding subsystems and values
    topSubsystems = T_FVA_ONOFF.uniSys(topIndex);
    heatmapData = [res2_keep(topIndex, comp1), res2_keep(topIndex, comp2)];

    % Generate the heatmap
    figure;
    heatmap(modelNames([comp1, comp2]), topSubsystems, heatmapData, ...
        'Colormap', custom_colormap, ...
        'ColorLimits', [min(heatmapData(:)), max(heatmapData(:))], ...
        'ColorbarVisible', 'on');

    title(['Top 41 Different Pathways: ', modelNames{comp1}, ' vs. ', modelNames{comp2}]);
    xlabel('Cancer Models');
    ylabel('Pathways');

    % Save the heatmap
    saveas(gcf, ['FVA_similarity_heatmap_', modelNames{comp1}, '_vs_', modelNames{comp2}, '.png']);
    saveas(gcf, ['FVA_similarity_heatmap_', modelNames{comp1}, '_vs_', modelNames{comp2}, '.pdf']);
end

%%
KICH_decrese = readtable('KICH_SFVA41.xlsx', 'PreserveVariableNames', true);
KIRC_decrese = readtable('KIRC_SFVA41.xlsx', 'PreserveVariableNames', true);
KIRP_decrese = readtable('KIRP_SFVA41.xlsx', 'PreserveVariableNames', true);
LIHC_decrese = readtable('LIHC_SFVA41.xlsx', 'PreserveVariableNames', true);

matSignalTKICH=KICH_decrese(1:30,2);
matSignalTKIRC=KIRC_decrese(1:30,2);
matSignalTKIRP=KIRP_decrese(1:30,2);
matSignalTLIHC=LIHC_decrese(1:30,2);

% matSignal=[matSignalT7,matSignalTH];
% matSignal=table2array(matSignal);
custom_colormap = redbluecmap(103); % Use redbluecmap function
unisubSys_TKICH= KICH_decrese(1:30,1);
unisubSys_TKICH=table2array(unisubSys_TKICH);
matSignalTKICH=table2array(matSignalTKICH);
%matSignalTH=table2array(matSignalTH);

%unisubSys_name_H= T_FVA_H_decrese.unisubSys(1:30);

name_7_RCC=strcat('Similarity_FVA_KICH');
%name_H_HCC=strcat('Similarity_FVA_all_subsystems_LIVER_CANCER');
modelNames_7_RCC = {'KICH\_High\_Low'};
%modelNames_H_HCC = {'LIHC\_ON\_OFF'};
modelNames_7_decrease = {'model_KICH'}; 
%modelNames_H_decrease = {'model_LIHC'};

modelNames_7_decrease2 = {'model_KICH',}; 
%modelNames_H_decrease2 = {'model_LIHC'};

feature astheightlimit 2000
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar


%%
%Create the clustergram
cgo_J_7 = clustergram(matSignalTKICH, ...
    'RowLabels', unisubSys_TKICH, ...
    'ColumnLabels', modelNames_7_decrease, ...
    'ColumnLabelsRotate', 360, ... % Adjusted for better readability
    'Cluster', 'all', ...
    'symmetric', false, ...
    'Colormap', altcolor);
% Add a title to the clustergram
addTitle(cgo_J_7, 'Similarity based on Flux Variability Analysis KICH tumor');

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo_J_7); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'FVA_similarity_heatmap_KICH_.svg');
saveas(cgf, 'FVA_similarity_heatmap_KICH_.pdf');

%%
% cgo_J_H= clustergram( matSignalTH(:,:),...
%     'RowLabels', unisubSys_name_H,...
%     'ColumnLabels', modelNames_H_decrease,...
%     'ColumnLabelsRotate',360, ...
%     'Cluster', 'all', ...
%     'symmetric','False',...
%     'Colormap', altcolor);
% 
% % Add a title to the clustergram
% addTitle(cgo_J_H, 'Similarity based on Flux Variability Analysis Liver Cancer');
% 
% % Wait for the clustergram to be created
% drawnow; % Ensure that the clustergram is rendered
% 
% % Get handle to the clustergram figure
% cgf = plot(cgo_J_H); % This should be a figure handle
% 
% % Add colorbar to the figure
% colorbar(cgf,'eastoutside');
% 
% % Save the figure with color scale
% saveas(cgf, 'FVA_similarity_heatmap_Honoff_.png');
% saveas(cgf, 'FVA_similarity_heatmap_Honoff_.pdf');

%% ???????????????????? ask Maria!!!!!
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_Pyruvate=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'ROS detoxification')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'ROS detoxification')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'ROS detoxification')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'ROS detoxification')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'ROS detoxification')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_ROS=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
% Table_sum_nucleotide=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_OXO=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_PPP=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Urea cycle')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Urea cycle')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Urea cycle')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Urea cycle')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Urea cycle')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_urea=table(rxns,rxnNames, formula,genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
% genes=consistent_model.grRules(ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis'),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_glycolysis=table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Nucleotide sugar metabolism')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_Nucleotide_sugar_metabolism =table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Glutathione metabolism')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_Glutathione_metabolism =table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Citric acid cycle')),:);
% maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
% rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
% rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
% genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
% formula= printRxnFormula(consistent_model,rxns);
% Table_sum_Citric_acid_cycle =table(rxns,rxnNames,formula, genes, mini, maxi);
% clear mini maxi rxns rxnNames genes formulas
% 
% %% FBA analysis
% 
% FBAsolution_control = optimizeCbModel(model_H_ctrl,'max','zero');
% FBAsolution_control.f % OF FBA result
% FBAsolution_control.x % FBA results of all reactions
% FBA_tableH = table(model_H_ctrl.subSystems,model_H_ctrl.rxns,FBAsolution_control.x);
% 
% FBAsolution_si1 = optimizeCbModel(model_H_sc1,'max','zero');
% FBA_table_sc1H = table(model_H_sc1.subSystems,model_H_sc1.rxns,FBAsolution_si1.x);
% 
% FBAsolution_si2 = optimizeCbModel(model_H_sc2,'max','zero');
% FBA_table_sc2H = table(model_H_sc2.subSystems,model_H_sc2.rxns,FBAsolution_si2.x);
% 
% FBA_table_final = table(consistent_model.subSystems,consistent_model.rxns);
% FBA_table_final.Control = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.Control((ismember(consistent_model.rxns,model_H_ctrl.rxns))) = FBAsolution_control.x;
% FBA_table_final.Control = cellstr(FBA_table_final.Control);
% 
% FBA_table_final.si1 = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.si1((ismember(consistent_model.rxns,model_H_sc1.rxns))) = FBAsolution_si1.x;
% FBA_table_final.si1 = cellstr(FBA_table_final.si1);
% 
% FBA_table_final.si2 = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.si2((ismember(consistent_model.rxns,model_H_sc2.rxns))) = FBAsolution_si2.x;
% FBA_table_final.si2 = cellstr(FBA_table_final.si2);
% 
% match=sum(strcmp(table2array(FBA_table_final(:,3:5)),'NA'),2)==3;
% FBA_tableH_final2=FBA_table_final(~match,:);
% %%
% 
% FBAsolution_control = optimizeCbModel(model_7_ctrl,'max','zero');
% FBAsolution_control.f % OF FBA result
% FBAsolution_control.x % FBA results of all reactions
% FBA_table = table(model_7_ctrl.subSystems,model_7_ctrl.rxns,FBAsolution_control.x);
% 
% FBAsolution_si1 = optimizeCbModel(model_7_sc1,'max','zero');
% FBA_table_sc1 = table(model_7_sc1.subSystems,model_7_sc1.rxns,FBAsolution_si1.x);
% 
% FBAsolution_si2 = optimizeCbModel(model_7_sc2,'max','zero');
% FBA_table_sc2 = table(model_7_sc2.subSystems,model_7_sc2.rxns,FBAsolution_si2.x);
% 
% FBA_table_final = table(consistent_model.subSystems,consistent_model.rxns);
% FBA_table_final.Control = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.Control((ismember(consistent_model.rxns,model_7_ctrl.rxns))) = FBAsolution_control.x;
% FBA_table_final.Control = cellstr(FBA_table_final.Control);
% 
% FBA_table_final.si1 = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.si1((ismember(consistent_model.rxns,model_7_sc1.rxns))) = FBAsolution_si1.x;
% FBA_table_final.si1 = cellstr(FBA_table_final.si1);
% 
% FBA_table_final.si2 = repmat("NA",numel(consistent_model.rxns),1);
% FBA_table_final.si2((ismember(consistent_model.rxns,model_7_sc2.rxns))) = FBAsolution_si2.x;
% FBA_table_final.si2 = cellstr(FBA_table_final.si2);
% 
% match=sum(strcmp(table2array(FBA_table_final(:,3:5)),'NA'),2)==3;
% FBA_table7_final2=FBA_table_final(~match,:);
% 
%% run the single gene deletion for 769-P
load models_KO_GLO1.mat
for i=1:numel(modelNames_kidney)
    name =  modelNames_kidney(i);
    model = models_Kidney_RPMIbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_biomass),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRateWT_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('essential', modelNames_kidney(i));
    save (name2{1},'TEssential');

end
save singlegenedeletionKidney_KO
%% run the single gene deletion for Huh7
load models_KO_GLO1.mat
for i=1:numel(modelNames_liver)
    name =  modelNames_liver(i);
    model = models_Liver_DMEMbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_ATP),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('essential', modelNames_liver(i));
    save (name2{1},'TEssential');

end
save singlegenedeletionLiver_KO
save analysis_script2onoffall