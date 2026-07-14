clear; close all; clc;
feature astheightlimit 2000
solverOK=1;

% Main script for metabolic modeling analysis

%%
% Load or generate metabolic network model
load ('driverModel_WithoutO2','colnamesH','colnames7','consistent_model','models_keep_7sample','models_keep_Hsample');

%load modelNames_Hwithouto2s.mat
%load modelNames_7withouto2s.mat

load GeneDrugRelations.mat

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

feature astheightlimit 2000

%% Jaccard Similarity Analysis
feature astheightlimit 2000
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar

%%
J = squareform(pdist(models_keep_Hsample','jaccard'));
%Jaccard similarity plots for sample models H

cgo_J = clustergram(1-J,...
    'RowLabels', colnamesH,...
    'ColumnLabels', colnamesH,...
    'ColumnLabelsRotate',45, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo_J,{'Model similarity based on Jaccard distance','models_keep'})

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo_J); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'ModelsimilaritybasedonJaccarddistance_H.png');

%%
J = squareform(pdist(models_keep_7sample','jaccard'));
%Jaccard similarity plots for sample models 7

cgo_J = clustergram(1-J,...
    'RowLabels', colnames7,...
    'ColumnLabels', colnames7,...
    'ColumnLabelsRotate',45, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo_J,{'Model similarity based on Jaccard distance','models_keep'})

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo_J); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'ModelsimilaritybasedonJaccarddistance_7.png');

%% Pathway analysis 
Pathways = table(unique(consistent_model.subSystems));
[pathways, ~, ub] = unique(consistent_model.subSystems);
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.consistent(ia) = T.path_counts(ib);
Pathways.Properties.VariableNames{1}='Pathways';

%% 
PathwaysH=Pathways;
% pathway information for the consensus models for H 
for i=1: numel(colnamesH)
   
[pathways, ~, ub] = unique(consistent_model.subSystems(models_keep_Hsample(:,i)~=0));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(PathwaysH.Pathways, T.pathways);
PathwaysH.Var2(ia) = T.path_counts(ib) ;
PathwaysH.Properties.VariableNames{2+i} = colnamesH{i};
end

Pathways7=Pathways;
% pathway information for the consensus models for 7 
for i=1: numel(colnames7)
[pathways, ~, ub] = unique(consistent_model.subSystems(models_keep_7sample(:,i)~=0));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways7.Pathways, T.pathways);
Pathways7.Var2(ia) = T.path_counts(ib) ;
Pathways7.Properties.VariableNames{2+i} = colnamesH{i};
end

%% plotting pathways
PathwayActivityH = PathwaysH;
for i=3:size(PathwayActivityH,2)
    PathwayActivityH(:,i) = array2table(table2array(PathwayActivityH(:,i))./table2array(PathwayActivityH(:,2)));
end

cgo = clustergram(table2array(PathwayActivityH(:,3:end)),...
    'RowLabels', PathwayActivityH.Pathways,...
    'ColumnLabels', PathwayActivityH.Properties.VariableNames(3:end),...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
h = plot(cgo); set(h,'TickLabelInterpreter','none');
colorbar(h)
title(h,'Pathway activity for all models H')

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'Pathwayactivityforallmodels_H.png');


PathwayActivity7 = Pathways7;
for i=3:size(PathwayActivity7,2)
    PathwayActivity7(:,i) = array2table(table2array(PathwayActivity7(:,i))./table2array(PathwayActivity7(:,2)));
end

cgo = clustergram(table2array(PathwayActivity7(:,3:end)),...
    'RowLabels', PathwayActivity7.Pathways,...
    'ColumnLabels', PathwayActivity7.Properties.VariableNames(3:end),...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
h = plot(cgo); set(h,'TickLabelInterpreter','none');
colorbar(h)
title(h,'Pathway activity for all models 7')

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'Pathwayactivityforallmodels_7.png');

%% Single Gene Deletion Analysis
geneList_keep=struct();
%% run the single gene deletion for 769-P

for i=1:numel(models_7_RPMIbiomass)
    name =  modelNames_7(i);
    model = models_7_RPMIbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_biomass),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRateWT_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('KO_GLO1_essential', modelNames_7(i));
    save (name2{1},'TEssential');

end

%% run the single gene deletion for Huh7
for i=1:numel(modelNames_H)
    name =  modelNames_H(i);
    model = models_H_DMEMbiomass{i}; %for each model
    model = changeObjective(model,'DM_atp_c_');

    [grRatio_ATP, grRateKO_ATP, grRateWT_ATP, hasEffect, ~, ~, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
   
   grRateWT_ATP=ones(numel(grRatio_ATP),1)*grRateWT_ATP;
   
    %Biomass
   
    model = changeObjective(model,'biomass_reaction');

     [grRatio_biomass, grRateKO_biomass, grRateWT_biomass, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(model,'FBA',[],0,1);
       grRateWT_biomass=ones(numel(grRatio_ATP),1)*grRateWT_biomass;

     TEssential= table(geneList,grRatio_biomass,grRateKO_biomass,grRatio_ATP,grRateWT_ATP,grRateKO_ATP)
    name2= strcat('KO_GLO1_essential', modelNames_H(i));
    save (name2{1},'TEssential');

end

%% Drug deletion
% define a list of drugs 
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
load GeneDrugRelations.mat
DrugList = unique(GeneDrugRelations.DrugName);
model_7_ctrl= changeObjective(model_7_ctrl,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_ctrl,'FBA',DrugList);
    
Drug_grRatio_biomass7(:,1)    = grRatio;
Drug_grRateKO_biomass7(:,1)   = grRateKO;
Drug_grRateWT_biomass7(:,1)   = grRateWT;
    
control_7_drug= DrugList(Drug_grRatio_biomass7(:,1) ==0);
save control_7_drug control_7_drug

%%
model_7_sc1= changeObjective(model_7_sc1,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc1,'FBA',DrugList);
    
Drug_grRatio_biomass7(:,2)    = grRatio;
Drug_grRateKO_biomass7(:,2)   = grRateKO;
Drug_grRateWT_biomass7(:,2)   = grRateWT;
    
Sc1_7_Drug = DrugList(Drug_grRatio_biomass(:,2) ==0);
save Sc1_7_Drug Sc1_7_Drug
%%
model_7_sc2= changeObjective(model_7_sc2,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc2,'FBA',DrugList);
    
Drug_grRatio_biomass7(:,3)    = grRatio;
Drug_grRateKO_biomass7(:,3)   = grRateKO;
Drug_grRateWT_biomass7(:,3)   = grRateWT;

    
Sc2_7_Drug = DrugList(Drug_grRatio_biomass(:,3) ==0);
save Sc2_7_Drug Sc2_7_Drug

%%
model_7_sc12= changeObjective(model_7_sc12,'biomass_reaction');   
    
[grRatio, grRateKO, grRateWT] = DrugDeletion(model_7_sc12,'FBA',DrugList);
    
Drug_grRatio_biomass7(:,4)    = grRatio;
Drug_grRateKO_biomass7(:,4)   = grRateKO;
Drug_grRateWT_biomass7(:,4)   = grRateWT;
    
Sc12_7_drug = DrugList(Drug_grRatio_biomass7(:,4) ==0);
save Sc12_7_drug Sc12_7_drug 

save workspace_singleDrugDeletion

%% Similarity Based on Flux Variability Analysis
%FVA 7
[minFlux_7contr, maxFlux_7contr] = fluxVariability(model_7_ctrl);
[minFlux_7sc1, maxFlux_7sc1] = fluxVariability(model_7_sc1);   
[minFlux_7sc12, maxFlux_7sc12] = fluxVariability(model_7_sc12); 
[minFlux_7sc2, maxFlux_7sc2] = fluxVariability(model_7_sc2);    

T_fluxes_7ctrl=table(model_7_ctrl.rxns,model_7_ctrl.rxnNames,minFlux_7contr,maxFlux_7contr);
T_fluxes_7sc1=table(model_7_sc1.rxns,model_7_sc1.rxnNames,minFlux_7sc1,maxFlux_7sc1);
T_fluxes_7sc2=table(model_7_sc2.rxns,model_7_sc2.rxnNames,maxFlux_7sc2);
T_fluxes_7sc12=table(model_7_sc12.rxns,model_7_sc12.rxnNames,maxFlux_7sc12);

[minFlux_Hcontr, maxFlux_Hcontr] = fluxVariability(model_H_ctrl);
[minFlux_Hsc1, maxFlux_Hsc1] = fluxVariability(model_H_sc1);
[minFlux_Hsc2, maxFlux_Hsc2] = fluxVariability(model_H_sc2);
[minFlux_Hsc12, maxFlux_Hsc12] = fluxVariability(model_H_sc12);

T_fluxes_Hctrl=table(model_H_ctrl.rxns,model_H_ctrl.rxnNames,minFlux_Hcontr,maxFlux_Hcontr);
T_fluxes_Hsc1=table(model_H_sc1.rxns,model_H_sc1.rxnNames,minFlux_Hsc1,maxFlux_Hsc1);
T_fluxes_Hsc2=table(model_H_sc2.rxns,model_H_sc2.rxnNames,maxFlux_Hsc2);
T_fluxes_Hsc12=table(model_H_sc12.rxns,model_H_sc12.rxnNames,maxFlux_Hsc12);

save analysis_part1

minFlux_7_keep=zeros(numel(consistent_model.rxns),4);
maxFlux_7_keep=zeros(numel(consistent_model.rxns),4);
minFlux_7_keep(ismember(consistent_model.rxns,model_7_ctrl.rxns),1)=minFlux_7contr;
maxFlux_7_keep(ismember(consistent_model.rxns,model_7_ctrl.rxns),1)=maxFlux_7contr;
minFlux_7_keep(ismember(consistent_model.rxns,model_7_sc1.rxns),2)=minFlux_7sc1;
maxFlux_7_keep(ismember(consistent_model.rxns,model_7_sc1.rxns),2)=maxFlux_7sc1;
minFlux_7_keep(ismember(consistent_model.rxns,model_7_sc2.rxns),3)=minFlux_7sc2;
maxFlux_7_keep(ismember(consistent_model.rxns,model_7_sc2.rxns),3)=maxFlux_7sc2;
minFlux_7_keep(ismember(consistent_model.rxns,model_7_sc12.rxns),3)=minFlux_7sc12;
maxFlux_7_keep(ismember(consistent_model.rxns,model_7_sc12.rxns),3)=maxFlux_7sc12;
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
T_FVA_7=table(uniSys,res2_keep);

save('T_FVA_7.mat','T_FVA_7')  % function form

%FVA huh7
minFlux_H_keep=zeros(numel(consistent_model.rxns),4);
maxFlux_H_keep=zeros(numel(consistent_model.rxns),4);
minFlux_H_keep(ismember(consistent_model.rxns,model_H_ctrl.rxns),1)=minFlux_Hcontr;
maxFlux_H_keep(ismember(consistent_model.rxns,model_H_ctrl.rxns),1)=maxFlux_Hcontr;
minFlux_H_keep(ismember(consistent_model.rxns,model_H_sc1.rxns),2)=minFlux_Hsc1;
maxFlux_H_keep(ismember(consistent_model.rxns,model_H_sc1.rxns),2)=maxFlux_Hsc1;
minFlux_H_keep(ismember(consistent_model.rxns,model_H_sc2.rxns),3)=minFlux_Hsc2;
maxFlux_H_keep(ismember(consistent_model.rxns,model_H_sc2.rxns),3)=maxFlux_Hsc2;
minFlux_H_keep(ismember(consistent_model.rxns,model_H_sc12.rxns),3)=minFlux_Hsc12;
maxFlux_H_keep(ismember(consistent_model.rxns,model_H_sc12.rxns),3)=maxFlux_Hsc12;
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
res2_keep=zeros(1,3);
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
T_FVA_H=table(uniSys,res2_keep);

save('T_FVA_H.mat','T_FVA_H')  
%% heatmap_FVAsimilarity

T_FVA_7_decrese = readtable('FVA_7_decrese_sc12.xlsx', 'PreserveVariableNames', true);
T_FVA_H_decrese = readtable('FVA_H_decrese_sc12.xlsx', 'PreserveVariableNames', true);

matSignalT7=T_FVA_7_decrese(1:30,2:4);
matSignalTH=T_FVA_H_decrese(1:30,2:4);
matSignal=[matSignalT7,matSignalTH];
matSignal=table2array(matSignal);
custom_colormap = redbluecmap(101); % Use redbluecmap function

matSignalT7=table2array(matSignalT7);
matSignalTH=table2array(matSignalTH);
unisubSys_name_7= T_FVA_7_decrese.unisubSys(1:30);
unisubSys_name_H= T_FVA_H_decrese.unisubSys(1:30);

name_7_RCC=strcat('Similarity_FVA_all_subsystems_RCC');
name_H_HCC=strcat('Similarity_FVA_all_subsystems_HCC');
modelNames_7_RCC = {'RCC\_ctrl\_sc1', 'RCC\_ctrl\_sc2', 'RCC\_ctrl\_sc12'};
modelNames_H_HCC = {'HCC\_ctrl\_sc1', 'HCC\_ctrl\_sc2', 'HCC\_ctrl\_sc12'};
modelNames_7_decrease = {'model_7_sc1', 'model_7_sc2', 'model_7_sc12'}; 
modelNames_H_decrease = {'model_H_sc1', 'model_H_sc2', 'model_H_sc12'};

feature astheightlimit 2000
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar


%%
% Create the clustergram
cgo_J_7 = clustergram(matSignalT7, ...
    'RowLabels', unisubSys_name_7, ...
    'ColumnLabels', modelNames_7_decrese, ...
    'ColumnLabelsRotate', 360, ... % Adjusted for better readability
    'Cluster', 'all', ...
    'symmetric', false, ...
    'Colormap', altcolor);
% Add a title to the clustergram
addTitle(cgo_J_7, 'Similarity based on Flux Variability Analysis RCC');

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo_J_7); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'FVA_similarity_heatmap_7_.png');

%%
cgo_J_H= clustergram( matSignalTH(:,:),...
    'RowLabels', unisubSys_name_H,...
    'ColumnLabels', modelNames_H_decrease,...
    'ColumnLabelsRotate',360, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);

% Add a title to the clustergram
addTitle(cgo_J_H, 'Similarity based on Flux Variability Analysis HCC');

% Wait for the clustergram to be created
drawnow; % Ensure that the clustergram is rendered

% Get handle to the clustergram figure
cgf = plot(cgo_J_H); % This should be a figure handle

% Add colorbar to the figure
colorbar(cgf,'eastoutside');

% Save the figure with color scale
saveas(cgf, 'FVA_similarity_heatmap_H_.png');

%% ???????????????????? ask Maria!!!!!
mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Pyruvate metabolism')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_Pyruvate=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'ROS detoxification')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'ROS detoxification')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'ROS detoxification')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'ROS detoxification')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'ROS detoxification')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_ROS=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Nucleotide metabolism')),:);
Table_sum_nucleotide=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Oxidative phosphorylation')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_OXO=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Pentose phosphate pathway')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_PPP=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Urea cycle')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Urea cycle')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Urea cycle')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Urea cycle')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Urea cycle')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_urea=table(rxns,rxnNames, formula,genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis')),:);
genes=consistent_model.grRules(ismember(consistent_model.subSystems,'Glycolysis/gluconeogenesis'),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_glycolysis=table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Nucleotide sugar metabolism')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Nucleotide sugar metabolism')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_Nucleotide_sugar_metabolism =table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Glutathione metabolism')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Glutathione metabolism')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_Glutathione_metabolism =table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

mini=minFlux_H_keep((ismember(consistentmodel.subSystems,'Citric acid cycle')),:);
maxi=maxFlux_H_keep((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
rxns=consistent_model.rxns((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
rxnNames=consistent_model.rxnNames((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
genes=consistent_model.grRules((ismember(consistent_model.subSystems,'Citric acid cycle')),:);
formula= printRxnFormula(consistent_model,rxns);
Table_sum_Citric_acid_cycle =table(rxns,rxnNames,formula, genes, mini, maxi);
clear mini maxi rxns rxnNames genes formulas

%% FBA analysis

FBAsolution_control = optimizeCbModel(model_H_ctrl,'max','zero');
FBAsolution_control.f % OF FBA result
FBAsolution_control.x % FBA results of all reactions
FBA_tableH = table(model_H_ctrl.subSystems,model_H_ctrl.rxns,FBAsolution_control.x);

FBAsolution_si1 = optimizeCbModel(model_H_sc1,'max','zero');
FBA_table_sc1H = table(model_H_sc1.subSystems,model_H_sc1.rxns,FBAsolution_si1.x);

FBAsolution_si2 = optimizeCbModel(model_H_sc2,'max','zero');
FBA_table_sc2H = table(model_H_sc2.subSystems,model_H_sc2.rxns,FBAsolution_si2.x);

FBA_table_final = table(consistent_model.subSystems,consistent_model.rxns);
FBA_table_final.Control = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.Control((ismember(consistent_model.rxns,model_H_ctrl.rxns))) = FBAsolution_control.x;
FBA_table_final.Control = cellstr(FBA_table_final.Control);

FBA_table_final.si1 = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.si1((ismember(consistent_model.rxns,model_H_sc1.rxns))) = FBAsolution_si1.x;
FBA_table_final.si1 = cellstr(FBA_table_final.si1);

FBA_table_final.si2 = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.si2((ismember(consistent_model.rxns,model_H_sc2.rxns))) = FBAsolution_si2.x;
FBA_table_final.si2 = cellstr(FBA_table_final.si2);

match=sum(strcmp(table2array(FBA_table_final(:,3:5)),'NA'),2)==3;
FBA_tableH_final2=FBA_table_final(~match,:);
%%

FBAsolution_control = optimizeCbModel(model_7_ctrl,'max','zero');
FBAsolution_control.f % OF FBA result
FBAsolution_control.x % FBA results of all reactions
FBA_table = table(model_H_ctrl.subSystems,model_7_ctrl.rxns,FBAsolution_control.x);

FBAsolution_si1 = optimizeCbModel(model_7_sc1,'max','zero');
FBA_table_sc1 = table(model_7_sc1.subSystems,model_7_sc1.rxns,FBAsolution_si1.x);

FBAsolution_si2 = optimizeCbModel(model_7_sc2,'max','zero');
FBA_table_sc2 = table(model_7_sc2.subSystems,model_7_sc2.rxns,FBAsolution_si2.x);

FBA_table_final = table(consistent_model.subSystems,consistent_model.rxns);
FBA_table_final.Control = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.Control((ismember(consistent_model.rxns,model_7_ctrl.rxns))) = FBAsolution_control.x;
FBA_table_final.Control = cellstr(FBA_table_final.Control);

FBA_table_final.si1 = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.si1((ismember(consistent_model.rxns,model_7_sc1.rxns))) = FBAsolution_si1.x;
FBA_table_final.si1 = cellstr(FBA_table_final.si1);

FBA_table_final.si2 = repmat("NA",numel(consistent_model.rxns),1);
FBA_table_final.si2((ismember(consistent_model.rxns,model_7_sc2.rxns))) = FBAsolution_si2.x;
FBA_table_final.si2 = cellstr(FBA_table_final.si2);

match=sum(strcmp(table2array(FBA_table_final(:,3:5)),'NA'),2)==3;
FBA_table7_final2=FBA_table_final(~match,:);

save analysis_script2