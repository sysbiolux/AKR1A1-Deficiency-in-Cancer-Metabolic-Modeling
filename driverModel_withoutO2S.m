%% loading model, medium, settings, dico
clear; close all; clc;
feature astheightlimit 2000
solverOK=1;

% model.rev = model.lb < 0; % the rev field will be a logical array
load consistent_model.mat
load dico2columns.mat
consistent_model.lb(ismember(consistent_model.rxns,'EX_o2s[e]'))=0;

% A = fastcc_4_rfastcormics(model, 1e-4, 0); % running FASTCC. 1e-4 needs to be used to avoid errors with small numbers.
already_mapped_tag = 0;
consensus_proportion = 0.9;
epsilon = 1e-4;

%% consensus model building H 
load discretizedH
rownames=geneSyms;

load medium_DMEM_H.mat
unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = consistent_model.rxns(ismember(consistent_model.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;
optional_settings.medium=medium_DMEM_H;
optional_settings.not_medium_constrained=[];
optional_settings.func = {'DM_atp_c_', 'biomass_reaction'}; % forced reactions

biomass_rxn='biomass_reaction';
% control model
[model_H_ctrl, ~] = fastcormics_RNAseq(consistent_model, discretized(:,[1,2,3,11]), rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% sc1 model
[model_H_sc1, ~] = fastcormics_RNAseq(consistent_model, discretized(:,4:6), rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);


% sc2 model
[model_H_sc2, ~] = fastcormics_RNAseq(consistent_model, discretized(:,7:10), rownames, dico, biomass_rxn,already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% sc12 model
[model_H_sc12, ~] = fastcormics_RNAseq(consistent_model, discretized(:,4:10), rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);


% check if biomass is still included
%% consensus model building 7
load discretized7
rownames=geneSyms;

load medium_RPMI_7.mat
unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = consistent_model.rxns(ismember(consistent_model.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;
optional_settings.medium=medium_RPMI_7;
optional_settings.not_medium_constrained=[];
optional_settings.func = {'DM_atp_c_', 'biomass_reaction'}; % forced reactions

% control model
[model_7_ctrl, ~] = fastcormics_RNAseq(consistent_model, discretized(:,[1,2,3,4]), rownames, dico,biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% sc1 model
[model_7_sc1, ~] = fastcormics_RNAseq(consistent_model, discretized(:,[5,6,7]), rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% sc2 model
[model_7_sc2, ~] = fastcormics_RNAseq(consistent_model, discretized(:,[8,9,10,11]), rownames, dico,biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);

% sc12 model
[model_7_sc12, ~] = fastcormics_RNAseq(consistent_model, discretized(:,5:11), rownames, dico, biomass_rxn,already_mapped_tag, consensus_proportion, epsilon, optional_settings);

%% sample specific model building H 
load discretizedH
rownames=geneSyms;

medium_DMEM_H = readtable("medium_DMEM_H.xlsx", 'ReadVariableNames', false);
medium_DMEM_H = table2array(medium_DMEM_H);
save medium_DMEM_H medium_DMEM_H
load medium_DMEM_H.mat

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = consistent_model.rxns(ismember(consistent_model.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;
optional_settings.medium=medium_DMEM_H;
optional_settings.not_medium_constrained=[];
optional_settings.func = {'DM_atp_c_', 'biomass_reaction'}; % forced reactions

models_keep_Hsample = zeros(numel(consistent_model.rxns),numel(colnames7));
% one model at a time
for i = 1:numel(colnamesH) %for each sample
    disp('H: '), 
    disp(i)
    [ContextModel, A_keep] = fastcormics_RNAseq(consistent_model, discretized(:,i), rownames, dico,biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
     name=strcat('Sample',colnamesH{i});
    save (name, 'ContextModel');
    models_keep_Hsample(A_keep,i)=1; % Only use this matrix for structure analysis
end


%% sample specific model building 7
load discretized7
rownames=geneSyms;

medium_RPMI_7 = readtable("medium_RPMI_7.xlsx", 'ReadVariableNames', false);
medium_RPMI_7 = table2array( medium_RPMI_7);
save medium_RPMI_7 medium_RPMI_7
load medium_RPMI_7.mat

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};

unpenalized = consistent_model.rxns(ismember(consistent_model.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;
optional_settings.medium=medium_RPMI_7;
optional_settings.not_medium_constrained=[];
models_keep_7sample = zeros(numel(consistent_model.rxns),numel(colnames7));

% one model at a time
for i = 1:numel(colnames7) %for each sample
    disp('7: '),
    disp(i)
    [ContextModel, A_keep] = fastcormics_RNAseq(consistent_model, discretized(:,i), rownames, dico,biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
    name=strcat('Sample',colnames7{i});
    save (name, 'ContextModel');
    models_keep_7sample(A_keep,i)=1; % Only use this matrix for structure analysis
end

save driverModel_WithoutO2

%save ("drivermodels.mat", "model_7_ctrl", "model_7_sc1", "model_7_sc12", "model_7_sc2", "model_H_ctrl" ,"model_H_sc1","model_H_sc12","model_H_sc2");
