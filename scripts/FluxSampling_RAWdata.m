%% FluxSampling_RAW.m
% Flux sampling analysis of AKR1A1-deficient metabolic models.
% This script performs flux sampling comparison analysis for:
%   - 769-P models
%   - HuH7 models
%% Required input files:
%   consistent_model.mat
%   SamplingResults_medium_1500_model_1.mat ... model_8.mat
%   samplingResults_medium_KO_GLO1_1500_model_1.mat ... model_8.mat
%
% This script compares feasible flux distributions between control and
% AKR1A1-deficient genome-scale metabolic models.
%
% The workflow performs:
%   1. Loading of previously sampled metabolic models.
%   2. Flux Balance Analysis (FBA).
%   3. Flux Variability Analysis (FVA).
%   4. Mapping of sampled fluxes to the reference Recon3D model.
%   5. Statistical comparison of flux distributions.
%   6. Identification of altered metabolic reactions and pathways.
%   7. Export of complete and pathway-specific result tables.
%
% Input:
%   consistent_model.mat
%       Reference consistent metabolic reconstruction.
%
%   SamplingResults_*.mat
%       Flux sampling results generated for each condition.
%
% Output:
%   Excel tables containing:
%       - reaction flux statistics
%       - fold-change values
%       - significance scores
%       - pathway-specific results
%
% Requirements:
%   MATLAB 
%   COBRA Toolbox
%   IBM CPLEX solver
%
% Associated study:
%   AKR1A1 deficiency metabolic modeling workflow.
%%

clear
close all
clc

feature astheightlimit 2000;
changeCobraSolver('ibm_cplex');


% Run flux sampling analysis for each cell line.
performAnalysisGSE1500(1,2:4,'769-P');
performAnalysisGSE1500(5,6:8,'Huh7');

%% Run TCGA cancer-type comparisons
% Model pairs:
%   KIRC: model 1 versus model 2
%   KIRP: model 3 versus model 4
%   KICH: model 5 versus model 6
%   LIHC: model 7 versus model 8

function performAnalysisGSE1500(ctrl,treatments,model_prefix_str)

%% Perform flux sampling comparison analysis
% Inputs:
%   ctrl             - control model index
%   treatments       - AKR1A1-deficient model indices
%   model_prefix_str - cell line identifier ('769-P' or 'Huh7')
%
% Steps:
%   1. Load sampled metabolic models.
%   2. Perform FBA and FVA.
%   3. Map results to the reference model.
%   4. Calculate flux statistics.
%   5. Export reaction-level tables.


%% Load reference model
% Load the consistent reference reconstruction used to align reactions
% across all sampled TCGA models.

load('consistent_model.mat')

model_orig = consistent_model;

cutoffFluxDiff = 100;

rxnsOfInterest = model_orig.rxns;

%% Standardize selected reaction identifiers
% Harmonize selected reaction IDs before mapping sampled-model reactions
% to the reference reconstruction.
% Standardize reaction identifiers to match the reference model.
model_orig.rxns(find(contains(model_orig.rxns,'CREATt4_2_r'))) = cellstr('CREATt4_2_r');
model_orig.rxns(find(contains(model_orig.rxns,'CREATt4_2'))) = cellstr('CREATt4_2');
model_orig.rxns(find(contains(model_orig.rxns,'ABUTt4_2'))) = cellstr('ABUTt4_2');
model_orig.rxns(find(contains(model_orig.rxns,'SRTNt6_2'))) = cellstr('SRTNt6_2');
model_orig.rxns(find(contains(model_orig.rxns,'GLYBt4_2'))) = cellstr('GLYBt4_2');
%%
%% Process control-versus-treatment model pairs
% The outer loops iterate through the selected model indices.
% The inner loop loads both the medium and KO_GLO1 sampling conditions.
%%
for i = 1:numel(ctrl)
    for j = 1:numel(treatments)
        treatmentNames = {'medium','KO_GLO1'};
        % Construct input filenames for the current model pair and condition.
        for k = 1:numel(treatmentNames)
            treatmentName = treatmentNames{k};
            % Load flux-sampling results for the standard medium condition.
            if treatmentName == "medium"
                disp (treatmentName)
                fileA_medium = ['SamplingResults_medium_1500_model_' num2str(ctrl(i)) '.mat'];
                fileB_medium = ['SamplingResults_medium_1500_model_' num2str(treatments(j)) '.mat'];
                % Load data
                load(fileA_medium, 'x');
                model1_medium = x.modelSampling;
                data1_medium = x.samples(:, 1:1500);
                load(fileB_medium, 'x');
                model2_medium = x.modelSampling;
                data2_medium = x.samples(:, 1:1500);
                disp(['... model & data loading done for control ' num2str(i) ' and treatment ' num2str(j) ' cellLine medium' model_prefix_str]);
            % Load flux-sampling results for the GLO1-knockout condition.
            elseif treatmentName == "KO_GLO1"
                disp (treatmentName)
                fileA_KO = ['samplingResults_medium_' treatmentName '_1500_model_' num2str(ctrl(i)) '.mat'];
                fileB_KO = ['samplingResults_medium_' treatmentName '_1500_model_' num2str(treatments(j)) '.mat'];
                % Load data
                load(fileA_KO, 'x');
                model1_KO = x.modelSampling;
                data1_KO = x.samples(:, 1:1500);
                load(fileB_KO, 'x');
                model2_KO = x.modelSampling;
                data2_KO = x.samples(:, 1:1500);
                
                disp(['... model & data loading done for control ' num2str(i) ' and treatment ' num2str(j) ' cellLine KO' model_prefix_str]);
                
                
            else
                disp("end loading")
            end
        end
        
        
        %% medium
        %% FBA and FVA: standard medium condition
        % Set the biomass reaction as the model objective for both models.
        % Perform analysis for model 1,2
        % Set objective for biomass reaction
        model1_medium.c = zeros(numel(model1_medium.c), 1);
        model1_medium.c(ismember(model1_medium.rxns, 'biomass_reaction')) = 1;
        
        model2_medium.c = zeros(numel(model2_medium.c), 1);
        model2_medium.c(ismember(model2_medium.rxns, 'biomass_reaction')) = 1;
        
        % Perform FBA
        FBAsolution_A_medium = optimizeCbModel(model1_medium, 'max', 'zero');
        x_A_medium = FBAsolution_A_medium.x;
        f_A_medium = FBAsolution_A_medium.f;
        f0_A_medium = FBAsolution_A_medium.f0;
        f1_A_medium = FBAsolution_A_medium.f1;
        f2_A_medium = FBAsolution_A_medium.f2;
        
        FBAsolution_B_medium = optimizeCbModel(model2_medium, 'max', 'zero');
        x_B_medium = FBAsolution_B_medium.x;
        f_B_medium = FBAsolution_B_medium.f;
        f0_B_medium = FBAsolution_B_medium.f0;
        f1_B_medium = FBAsolution_B_medium.f1;
        f2_B_medium = FBAsolution_B_medium.f2;
        
        fba_medium=[];
        fba_medium=[fba_medium; f_A_medium, f_B_medium, f0_A_medium, f0_B_medium, f1_A_medium,  f1_B_medium, f2_A_medium, f2_B_medium];
        file_fba_medium=[model_prefix_str '_fba_medium' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
        if exist(file_fba_medium, 'file'), delete(file_fba_medium); end
        writematrix(fba_medium,file_fba_medium)
        disp(file_fba_medium)
        
        
         % Perform FVA
         [minFlux_A_medium, maxFlux_A_medium] = fluxVariability(model1_medium, 90); % Adjust the number of samples
         [minFlux_B_medium, maxFlux_B_medium] = fluxVariability(model2_medium, 90); % Adjust the number of samples
        
        %% KO
        %% FBA and FVA: GLO1-knockout condition
        % Set the biomass reaction as the model objective for both models.
        % Perform analysis for model 1,2
        % Set objective for biomass reaction
        model1_KO.c = zeros(numel(model1_KO.c), 1);
        model1_KO.c(ismember(model1_KO.rxns, 'biomass_reaction')) = 1;
        
        model2_KO.c = zeros(numel(model2_KO.c), 1);
        model2_KO.c(ismember(model2_KO.rxns, 'biomass_reaction')) = 1;
        
        % Perform FBA
        FBAsolution_A_KO = optimizeCbModel(model1_KO, 'max', 'zero'); %ctrl with KO GLO1
        x_A_KO = FBAsolution_A_KO.x;
        f_A_KO = FBAsolution_A_KO.f;
        f0_A_KO = FBAsolution_A_KO.f0;
        f1_A_KO = FBAsolution_A_KO.f1;
        f2_A_KO = FBAsolution_A_KO.f2;
        
        FBAsolution_B_KO = optimizeCbModel(model2_KO, 'max', 'zero');
        x_B_KO = FBAsolution_B_KO.x;
        f_B_KO = FBAsolution_B_KO.f;
        f0_B_KO = FBAsolution_B_KO.f0;
        f1_B_KO = FBAsolution_B_KO.f1;
        f2_B_KO = FBAsolution_B_KO.f2;
        
        fba_KO=[];
        fba_KO=[fba_medium; f_A_KO, f_B_KO, f0_A_KO, f0_B_KO, f1_A_KO,  f1_B_KO, f2_A_KO, f2_B_KO];
        file_fba_KO=[model_prefix_str '_fba_KO' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
        if exist(file_fba_KO, 'file'), delete(file_fba_KO); end
        writematrix( fba_KO,file_fba_KO)
        
         % Perform FVA
         [minFlux_A_KO, maxFlux_A_KO] = fluxVariability(model1_KO, 90); % Adjust the number of samples
         [minFlux_B_KO, maxFlux_B_KO] = fluxVariability(model2_KO, 90); % Adjust the number of samples

    %% mapping to original model
    %% Map results to the reference model
    % Align sampled-model reactions with the reference reaction list.
    %
    % Model A:
    %   control/reference model
    %
    % Model B:
    %   comparison model
    %A
    [C1m,IA1m,IB1m] = intersect(model_orig.rxns,model1_medium.rxns,'stable'); %medium 1
    
    resA_sampling_medium=zeros(numel(model_orig.rxns),1500);
    resA_sampling_medium(IA1m,:)=data1_medium; %data 1 medium
    
    resA_fba_m=zeros(numel(model_orig.rxns), 1);
    resA_fba_m(IA1m)=x_A_medium;
    
    resA_minflux_m=zeros(numel(model_orig.rxns), 1);
    resA_minflux_m(IA1m)=minFlux_A_medium;
     
    resA_maxflux_m=zeros(numel(model_orig.rxns), 1);
    resA_maxflux_m(IA1m)=maxFlux_A_medium;
%     
    resOrig_subSys_m=cell(numel(model_orig.rxns),1);
    resOrig_subSys_m=model_orig.subSystems
    
    resOrig_rxns_m=cell(numel(model_orig.rxns),1);
    resOrig_rxns_m=model_orig.rxns
    
    resOrig_rxnNames_m=cell(numel(model_orig.rxns),1);
    resOrig_rxnNames_m=model_orig.rxnNames;
    
    resOrig_formula_m=cell(numel(model_orig.rxns),1);
    resOrig_formula_m=printRxnFormula(model_orig);
    
    %% B
    [C2m,IA2m,IB2m] = intersect(model_orig.rxns,model2_medium.rxns,'stable'); %medium 2
    [C2ko,IA2ko,IB2ko] = intersect(model_orig.rxns,model2_KO.rxns,'stable'); %KO 2
    
    resB_sampling_medium=zeros(numel(model_orig.rxns),1500);
    resB_sampling_medium(IA2m,:)=data2_medium;%data2
    
    resB_sampling_KO=zeros(numel(model_orig.rxns),1500);
    resB_sampling_KO(IA2ko,:)=data2_KO;%data 2
    
    resB_fba_m=zeros(numel(model_orig.rxns), 1);
    resB_fba_m(IA2m)=x_B_medium;
    
    resB_minflux_m=zeros(numel(model_orig.rxns), 1);
    resB_minflux_m(IA2m)=minFlux_B_medium;
     
    resB_maxflux_m=zeros(numel(model_orig.rxns), 1);
    resB_maxflux_m(IA2m)=maxFlux_B_medium;
    
    resB_fba_ko=zeros(numel(model_orig.rxns), 1);
    resB_fba_ko(IA2ko)= x_B_KO;
    
    resB_minflux_ko=zeros(numel(model_orig.rxns), 1);
    resB_minflux_ko(IA2ko)=minFlux_B_KO;
     
    resB_maxflux_ko=zeros(numel(model_orig.rxns), 1);
    resB_maxflux_ko(IA2ko)=maxFlux_B_KO;
    
    resOrig_subSys_ko=cell(numel(model_orig.rxns),1);
    resOrig_subSys_ko=model_orig.subSystems;
    
    resOrig_rxns_ko=cell(numel(model_orig.rxns),1);
    resOrig_rxns_ko=model_orig.rxns;
    
    resOrig_rxnNames_ko=cell(numel(model_orig.rxns),1);
    resOrig_rxnNames_ko=model_orig.rxnNames;
    
    resOrig_formula_ko=cell(numel(model_orig.rxns),1);
    resOrig_formula_ko=printRxnFormula(model_orig);
    
    %mean(resA_sampling_medium(find(ismember(model_orig.rxns,'PYK')),:));
    mean(resB_sampling_medium(find(ismember(model_orig.rxns,'PYK')),:));
    %mean(resA_sampling_KO(find(ismember(model_orig.rxns,'PYK')),:));
    mean(resB_sampling_KO(find(ismember(model_orig.rxns,'PYK')),:));
    
    %% statistical test medium
    stats_medium=[];
    direction_change_m = zeros(numel(model_orig.rxns), 1);
    for counter=1:size(resA_sampling_medium,1)
        A=resA_sampling_medium(counter,:); % control
        B=resB_sampling_medium(counter,:); % treatment AK1R1 deficiency model
        P=ranksum(A,B);
        % Check conditions and perform calculations
        if (A >= 0 & B >= 0) | (A <= 0 & B <= 0)
            % Mark a change in direction (difference is zero)
            direction_change_m(counter) = 0;
            % Take the absolute values of A and B
            A = abs(A);
            B = abs(B);
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
        elseif A <= 0 & B >= 0
            % Mark a change in direction
            direction_change_m(counter) = 1;
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
            snr =abs(snr);
        elseif A >= 0 & B <= 0
            % Mark a change in direction
            direction_change_m(counter) = -1;
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
            snr = abs(snr);
        end
        
        % Check if the mean of A and B is zero
        if A == 0
            % If mean(B) is zero, add 0.1 to the sum
            A1 = A + 0.1;
        else
            % If mean(B) is not zero, just calculate the sum
            A1 = A;
        end
        % Check if the mean of B is zero
        if B == 0
            % If mean(B) is zero, add 0.1 to the sum
            B1 = B + 0.1;
        else
            % If mean(B) is not zero, just calculate the sum
            B1 = B;
        end
        %stats=[stats; mean(A), mean(B), mean(B)-mean(A), mean(B)/mean(A), P -log10(P)];
        log2FC=log2(abs(mean(B1)/mean(A1)));
        stats_medium=[stats_medium; mean(A), mean(B), std(A), std(B), mean(B)-mean(A), log2FC, snr, P, -log2(P)]; % Fold change with abs (module) to not have i numbers after log2
        
    end
    
    stats_medium(10:20,:)
    
    figure
    histogram(stats_medium(:,9))
    %title('P values (-log10)')
    title([model_prefix_str 'P values (-log2)'])
    
    figure
    %hist(log10(stats(:,4)))
    histogram(stats_medium(:,6))
    %title('log10 foldchange (mean(B)/mean(A))')
    title([model_prefix_str 'log2 foldchange (mean(B)/mean(A))'])
    
    figure
    % plot(log10(stats(:,4)),stats(:,6),'*')
    % title('vulcano: log10 foldchange vs -log10(P)')
    plot(stats_medium(:,6),stats_medium(:,9),'*')
    title([model_prefix_str 'vulcano: log2 foldchange vs -log2(P)'])
    
    statst_medium=array2table(stats_medium,'RowNames',model_orig.rxns,'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    % Create a table with columns A, B, and direction_change
    statst_medium.direction_change=direction_change_m;
    
    %% statistical test KO
    stats_KO=[];
    direction_change_ko = zeros(numel(model_orig.rxns), 1);
    for counter=1:size(resA_sampling_medium,1)
        A=resA_sampling_medium(counter,:);
        B=resB_sampling_KO(counter,:);
        P=ranksum(A,B);
        % Check conditions and perform calculations
        if (A >= 0 & B >= 0) | (A <= 0 & B <= 0)
            % Mark a change in direction (difference is zero)
            direction_change_ko(counter) = 0;
            % Take the absolute values of A and B
            A = abs(A);
            B = abs(B);
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
        elseif A <= 0 & B >= 0
            % Mark a change in direction
            direction_change_ko(counter) = 1;
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
            snr =abs(snr);
        elseif A >= 0 & B <= 0
            % Mark a change in direction
            direction_change_ko(counter) = -1;
            snr =(mean(B) - mean(A)) / (std(A) + std(B)); % signal to noise ratio
            snr = abs(snr);
        end
        % Check if the mean of A and B is zero
        if A == 0
            % If mean(B) is zero, add 0.1 to the sum
            A1 = A + 0.1;
        else
            % If mean(B) is not zero, just calculate the sum
            A1 = A;
        end
        % Check if the mean of B is zero
        if B == 0
            % If mean(B) is zero, add 0.1 to the sum
            B1 = B + 0.1;
        else
            % If mean(B) is not zero, just calculate the sum
            B1 = B;
        end
        %stats=[stats; mean(A), mean(B), mean(B)-mean(A), mean(B)/mean(A), P -log10(P)];
        %log2FC=log2(B1)-log2(A1);
        log2FC=log2(abs(mean(B1)/mean(A1)));
        stats_KO=[stats_KO; mean(A), mean(B), std(A), std(B), mean(B)-mean(A), log2FC, snr, P, -log2(P)];
    end
    
    
    stats_KO(1:10,:)
    
    
    figure
    histogram(stats_KO(:,9))
    %title('P values (-log10)')
    title([model_prefix_str 'P values (-log2)'])
    
    figure
    %hist(log10(stats(:,4)))
    histogram(stats_KO(:,6))
    %title('log10 foldchange (mean(B)/mean(A))')
    title([model_prefix_str 'log2 foldchange (mean(B)/mean(A))'])
    
    figure
    % plot(log10(stats(:,4)),stats(:,6),'*')
    % title('vulcano: log10 foldchange vs -log10(P)')
    plot(stats_KO(:,6),stats_KO(:,9),'*')
    title([model_prefix_str 'vulcano: log2 foldchange vs -log2(P)'])
    
    statst_KO=array2table(stats_KO,'RowNames',model_orig.rxns,'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    statst_KO.direction_change = direction_change_ko;
    
    %% writetable medium
    cutoffFluxDiff=0;
    %stats=setdiff(stats,cofactors);
    statst_medium.rxnNames=model_orig.rxnNames;
    statst_medium.subSystems=model_orig.subSystems;
    statst_medium.formula=resOrig_formula_m;
    %statst_medium.rxns=model_orig.rxns;
    file_medium=[model_prefix_str '_FluxSamplingStatst_all_medium_5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file_medium, 'file'), delete(file_medium); end
    writetable(statst_medium,file_medium,'WriteRowNames',true)
    
    rxnsOfInterest=model_orig.rxns;
    rxnsOfInterest=find(ismember(model_orig.rxns,rxnsOfInterest));
    
    %%'direction_change','rxnNames','subSystems','formula'
    
    dn=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    dn=intersect(dn,rxnsOfInterest);
    model_origformula=printRxnFormula(model_orig);
    dnt=array2table(stats_medium(dn,:),'RowNames',model_orig.rxns(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    dnt.direction_change = direction_change_m(dn);
    dnt.rxnNames=model_orig.rxnNames(dn);
    dnt.subSystems=model_orig.subSystems(dn);
    dnt.formula=model_origformula(dn);
    dnt=sortrows(dnt,6,'descend');
    
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    up=intersect(up,rxnsOfInterest);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.direction_change = direction_change_m(up);
    upt.rxnNames=model_orig.rxnNames(up);
    upt.subSystems=model_orig.subSystems(up);
    upt.formula=model_origformula(up);
    upt=sortrows(upt,6,'descend');
    
    file=[model_prefix_str '_FluxSamplingStatsallrank_medium_5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true,'Sheet','Up')
    writetable(dnt,file,'WriteRowNames',true,'Sheet','Down')
    
    
    %% writetable KO
    cutoffFluxDiff=0;
    %stats=setdiff(stats,cofactors);
    statst_KO.metNames=model_orig.rxnNames;
    statst_KO.subSystems=model_orig.subSystems;
    statst_KO.formula=printRxnFormula(model_orig);
    statst_KO.rxns=model_orig.rxns;
    
    file_KO=[model_prefix_str '_FluxSamplingStatst_all_KO_5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file_KO, 'file'), delete(file_KO); end
    writetable(statst_KO,file_KO,'WriteRowNames',true)
    
    dnko=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    dnko=intersect(dnko,rxnsOfInterest);
    dntko=array2table(stats_KO(dnko,:),'RowNames',model_orig.rxns(dnko),'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    dntko.direction_change = direction_change_ko(dnko);
    dntko.rxnNames=model_orig.rxnNames(dnko);
    dntko.subSystems=model_orig.subSystems(dnko);
    dntko.rxns=model_orig.rxns(dnko);
    dntko.formula=model_origformula(dnko);
    dntko=sortrows(dntko,6,'descend');
    
    upko=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    upko=intersect(upko,rxnsOfInterest);
    uptko=array2table(stats_KO(upko,:),'RowNames',model_orig.rxns(upko),'VariableNames',{'meanA','meanB','stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    uptko.direction_change = direction_change_ko(upko);
    uptko.rxnNames=model_orig.rxnNames(upko);
    uptko.subSystems=model_orig.subSystems(upko);
    uptko.rxns=model_orig.rxns(upko);
    uptko.formula=model_origformula(upko);
    uptko=sortrows(uptko,6,'descend');
    
    file=[model_prefix_str '_FluxSamplingStatsallrank_KO_5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(uptko,file,'WriteRowNames',true,'Sheet','Up')
    writetable(dntko,file,'WriteRowNames',true,'Sheet','Down')
    
    %% writetable by selected pathways
    %Glycolysis/gluconeogenesis
    cutoffFluxDiff=0;
    pathway='Glycolysis/gluconeogenesis';
    model_origformula=printRxnFormula(model_orig);
    
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList);
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    direction_change_up=direction_change_m(up);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    direction_change_up2=direction_change_m(up2);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    upt.formula=model_origformula(up);
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Glyco_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    direction_change_up=direction_change_ko(up);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    direction_change_up2=direction_change_ko(up2);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula;
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Glyco_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    %%
    %Pyruvate metabolism
    cutoffFluxDiff=0;
    pathway='Pyruvate metabolism';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList);
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Pyr_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Pyr_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    %%
    %Pentose phosphate pathway
    cutoffFluxDiff=0;
    pathway='Pentose phosphate pathway';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList)
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    
    file=[model_prefix_str '_fluxSamplingStats_a_PPP_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_PPP_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    %%
    %Glutathione metabolism
    cutoffFluxDiff=0;
    pathway='Glutathione metabolism';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList)
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    file=[model_prefix_str '_fluxSamplingStats_a_Glutath_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Glutath_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    %%
    %Cholesterol metabolism
    cutoffFluxDiff=0;
    pathway='Cholesterol metabolism';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList)
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    
    file=[model_prefix_str '_fluxSamplingStats_a_Chol_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Chol_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    %%
    %ROS detoxification metabolism
    cutoffFluxDiff=0;
    pathway='ROS detoxification';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList);
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    
    file=[model_prefix_str '_fluxSamplingStats_a_ROS_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_ROS_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    %%
    %TCA pathway
    cutoffFluxDiff=0;
    pathway='Citric acid cycle';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList)
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    
    file=[model_prefix_str '_fluxSamplingStats_a_TCA_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_TCA_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    %%
    %Nucleotide metabolism
    cutoffFluxDiff=0;
    pathway='Nucleotide interconversion';
    rxnList=find((ismember(model_orig.subSystems,pathway)));
    rxns=model_orig.rxns(rxnList)
    up=find((stats_medium(:,6)>0).*(stats_medium(:,9)>0).*(stats_medium(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_medium(:,6)<0).*(stats_medium(:,9)>0).*(stats_medium(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_medium(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    
    
    file=[model_prefix_str '_fluxSamplingStats_a_Nucleot_medium5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    
    up=find((stats_KO(:,6)>0).*(stats_KO(:,9)>0).*(stats_KO(:,5)>cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up=intersect(up,rxnList);
    up2=find((stats_KO(:,6)<0).*(stats_KO(:,9)>0).*(stats_KO(:,5)<-cutoffFluxDiff));
    %up=setdiff(up,cofactors);
    up2=intersect(up2,rxnList);
    up=union(up, up2);
    upt=array2table(stats_KO(up,:),'RowNames',model_orig.rxns(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','log2FC','SNR','pValue','-log2_p'});
    upt.rxnNames=model_orig.rxnNames(up);
    modelformula=  model_origformula(up);
    upt.formula=modelformula
    upt=sortrows(upt,7,'descend') %by SNR
    file=[model_prefix_str '_fluxSamplingStats_a_Nucleot_KO5' num2str(ctrl(i)) '_vs_' num2str(treatments(j)) '_1500.xlsx'];
    if exist(file, 'file'), delete(file); end
    writetable(upt,file,'WriteRowNames',true)
    end
end
end
