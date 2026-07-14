%% flux_sum_analysis.m
% Unified legacy Flux Sum analysis for:
%   1. Historical standard 1500-labelled workflow
%   2. TCGA ON/OFF 1000-sample workflow
%
%
% MATLAB compatibility: R2019

clearvars -except solverOK
clc
close all

%% Select one legacy mode
datasetType = 'standard1500';
% datasetType = 'onoff1000';

%% Relative folders
referenceFolder = '../data/reference';
standardSamplingFolder = '../data/sampling/standard';
onoffSamplingFolder = '../data/sampling/onoff';
standardOutputFolder = '../data/flux_sum/standard';
onoffOutputFolder = '../data/flux_sum/onoff';

referenceModelFile = fullfile(referenceFolder,'consistent_model.mat');
cofactorFile = fullfile(referenceFolder,'metsCofactors.txt');

load(referenceModelFile,'consistent_model');
model_orig = consistent_model;

epsilon = 1e-4; %#ok<NASGU>
metsOfInterest = model_orig.mets;
metsOfInterest = find(ismember(model_orig.mets,metsOfInterest)); %#ok<NASGU>

%% Select comparisons and folders
if strcmp(datasetType,'standard1500')
    comparisonA = [1 1 1 5 5 5];
    comparisonB = [3 2 4 7 6 8];
    comparisonPrefix = {'769-P','769-P','769-P','Huh7','Huh7','Huh7'};
    inputFolder = standardSamplingFolder;
    outputFolder = standardOutputFolder;
    nSamplesUsed = 500;
elseif strcmp(datasetType,'onoff1000')
    comparisonA = [1 3 5 7];
    comparisonB = [2 4 6 8];
    comparisonPrefix = {'KIRC-ONOFF','KIRP-ONOFF','KICH-ONOFF','LIHC-ONOFF'};
    inputFolder = onoffSamplingFolder;
    outputFolder = onoffOutputFolder;
    nSamplesUsed = 1000;
else
    error('datasetType must be ''standard1500'' or ''onoff1000''.');
end

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% Run all comparisons
for comparisonIndex = 1:numel(comparisonA)
    nrA = comparisonA(comparisonIndex);
    nrB = comparisonB(comparisonIndex);
    model_prefix_str = comparisonPrefix{comparisonIndex};

    if strcmp(datasetType,'standard1500')
        treatments = {'medium','KO_GLO1'};
        for treatmentIndex = 1:length(treatments)
            treatment = treatments{treatmentIndex};
            fileA = fullfile(inputFolder,['SamplingResults_medium_1500_model_' num2str(nrA) '.mat']);
            fileB = fullfile(inputFolder,['SamplingResults_medium_1500_model_' num2str(nrB) '.mat']);
            runFluxSumPairLegacy(model_orig,cofactorFile,fileA,fileB,nSamplesUsed, ...
                model_prefix_str,nrA,nrB,treatment,datasetType,outputFolder);
        end
    elseif strcmp(datasetType,'onoff1000')
        fileA = fullfile(inputFolder,['samplingResults_mediumonoff_1000_model_' num2str(nrA) '.mat']);
        fileB = fullfile(inputFolder,['samplingResults_mediumonoff_1000_model_' num2str(nrB) '.mat']);
        runFluxSumPairLegacy(model_orig,cofactorFile,fileA,fileB,nSamplesUsed, ...
            model_prefix_str,nrA,nrB,'medium',datasetType,outputFolder);

        fileA = fullfile(inputFolder,['samplingResults_medium_KO_GLO1onoff_1000_model_' num2str(nrA) '.mat']);
        fileB = fullfile(inputFolder,['samplingResults_medium_KO_GLO1onoff_1000_model_' num2str(nrB) '.mat']);
        runFluxSumPairLegacy(model_orig,cofactorFile,fileA,fileB,nSamplesUsed, ...
            model_prefix_str,nrA,nrB,'KO_GLO1',datasetType,outputFolder);
    end
end


function runFluxSumPairLegacy(model_orig,cofactorFile,fileA,fileB,nSamplesUsed, ...
    model_prefix_str,nrA,nrB,treatment,datasetType,outputFolder)

close all

disp(['Processing: ' model_prefix_str ' - ' treatment]);
disp(['Input A: ' fileA]);
disp(['Input B: ' fileB]);

load(fileA,'x');
model1 = x.modelSampling;
data1 = x.samples(:,1:nSamplesUsed);
size(data1)

load(fileB,'x');
model2 = x.modelSampling;
data2 = x.samples(:,1:nSamplesUsed);
size(data2)

disp('... model & data loading done ...')

res1 = [];
for counter = 1:size(data1,2)
    v = data1(:,counter);
    temp = repmat(v',size(model1.S,1),1);
    fluxes = model1.S.*temp;
    fluxSumP = full(sum((fluxes>0).*fluxes,2));
    fluxSumN = full(sum((fluxes<0).*fluxes,2)); %#ok<NASGU>
    res1 = [res1,fluxSumP]; %#ok<AGROW>
end

disp('... fluxSum A calculated ...')

res2 = [];
for counter = 1:size(data2,2)
    v = data2(:,counter);
    temp = repmat(v',size(model2.S,1),1);
    fluxes = model2.S.*temp;
    fluxSumP = full(sum((fluxes>0).*fluxes,2));
    fluxSumN = full(sum((fluxes<0).*fluxes,2)); %#ok<NASGU>
    res2 = [res2,fluxSumP]; %#ok<AGROW>
end

disp('... fluxSum B calculated ...')

figure
boxplot(res1(1:30,:)','Labels',model1.mets(1:30))
set(gca,'FontSize',10,'XTickLabelRotation',45)

if strcmp(datasetType,'onoff1000')
    figure
    boxplot(res2(1:30,:)','Labels',model2.mets(1:30))
    set(gca,'FontSize',10,'XTickLabelRotation',45)
end

mean(res1(find(ismember(model1.mets,'glc_D[c]')),:))
mean(res2(find(ismember(model2.mets,'glc_D[c]')),:))

[C,IA,IB] = intersect(model_orig.mets,model1.mets,'stable'); %#ok<ASGLU>
resA = zeros(numel(model_orig.mets),nSamplesUsed);
resA(IA,:) = res1;

[C,IA,IB] = intersect(model_orig.mets,model2.mets,'stable'); %#ok<ASGLU>
resB = zeros(numel(model_orig.mets),nSamplesUsed);
resB(IA,:) = res2;

mean(resA(find(ismember(model_orig.mets,'glc_D[c]')),:))
mean(resB(find(ismember(model_orig.mets,'glc_D[c]')),:))

stats = [];
for counter = 1:size(resA,1)
    A = resA(counter,:);
    B = resB(counter,:);
    P = ranksum(A,B);

    if A == 0
        A1 = A + 0.1;
    else
        A1 = A;
    end

    if B == 0
        B1 = B + 0.1;
    else
        B1 = B;
    end

    snr = (mean(B1)-mean(A1))/(std(A)+std(B))
    log2FC = log2(abs(mean(B1)/mean(A1)));

    stats = [stats; mean(A),mean(B),std(A),std(B), ...
        mean(B)-mean(A),log2FC,snr,P,-log2(P)]; %#ok<AGROW>
end

stats(1:10,:)

figure
histogram(stats(:,9))
title([model_prefix_str 'P values (-log2)'])

figure
histogram(stats(:,6))
title([model_prefix_str 'log2 foldchange (mean(B)/mean(A))'])

figure
plot(stats(:,6),stats(:,9),'*')
title([model_prefix_str 'vulcano: log2 foldchange vs -log2(P)'])

metsOfInterest = model_orig.mets;
metsOfInterest = find(ismember(model_orig.mets,metsOfInterest)); %#ok<NASGU>

str = fileread(cofactorFile);
cofactorNames = regexp(str,'\r\n|\r|\n','split')'
cofactors = find(ismember(model_orig.metNames,cofactorNames))

statst = array2table(stats,'RowNames',model_orig.mets, ...
    'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
statst.metNames = model_orig.metNames;

file = buildFluxSumOutputNameLegacy(outputFolder,model_prefix_str,'all', ...
    treatment,nrA,nrB,datasetType);
deleteIfPresentLegacy(file);
writetable(statst,file,'WriteRowNames',true);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'Glycolysis/gluconeogenesis','Glycolysis',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'Pentose phosphate pathway','PPP',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'Pyruvate metabolism','Pyruvate',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'Glutathione metabolism','Glutathione',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'ROS detoxification','ROS_detoxification',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

writeFluxSumPathwayLegacy(stats,model_orig,cofactors, ...
    'Citric acid cycle','Citric_acid_cycle',model_prefix_str,treatment, ...
    nrA,nrB,datasetType,outputFolder);

cutoffFluxDiff = 100;
metsOfInterest = model_orig.mets;
metsOfInterest = find(ismember(model_orig.mets,metsOfInterest));

up = find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
up = intersect(up,metsOfInterest);
upt = array2table(stats(up,:),'RowNames',model_orig.mets(up), ...
    'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
upt.metNames = model_orig.metNames(up);
upt = sortrows(upt,7,'descend');

dn = find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
dn = intersect(dn,metsOfInterest);
dnt = array2table(stats(dn,:),'RowNames',model_orig.mets(dn), ...
    'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
dnt.metNames = model_orig.metNames(dn);
dnt = sortrows(dnt,9,'descend');

file = buildFluxSumOutputNameLegacy(outputFolder,model_prefix_str,'allrank', ...
    treatment,nrA,nrB,datasetType);
deleteIfPresentLegacy(file);
writetable(upt,file,'WriteRowNames',true,'Sheet','Up');
writetable(dnt,file,'WriteRowNames',true,'Sheet','Down');

end


function writeFluxSumPathwayLegacy(stats,model_orig,cofactors,pathway, ...
    pathwayLabel,model_prefix_str,treatment,nrA,nrB,datasetType,outputFolder)

metsOfInterest = {};
metsOfInterest = model_orig.mets; %#ok<NASGU>
metsOfInterest = find(ismember(model_orig.mets,model_orig.mets)); %#ok<NASGU>

metList = findMetsFromRxns(model_orig, ...
    model_orig.rxns(ismember(model_orig.subSystems,pathway)));
metListidx = find(ismember(model_orig.mets,metList));

cutoffFluxDiff = 0;
up = find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
up = setdiff(up,cofactors);
up = intersect(up,metListidx);

up2 = find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
up2 = setdiff(up2,cofactors);
up2 = intersect(up2,metListidx);

up = union(up,up2);

upt = array2table(stats(up,:),'RowNames',model_orig.mets(up), ...
    'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
upt.metNames = model_orig.metNames(up);
upt = sortrows(upt,7,'descend');

file = buildFluxSumOutputNameLegacy(outputFolder,model_prefix_str,pathwayLabel, ...
    treatment,nrA,nrB,datasetType);
deleteIfPresentLegacy(file);
writetable(upt,file,'WriteRowNames',true);

end


function file = buildFluxSumOutputNameLegacy(outputFolder,modelPrefix,section, ...
    treatment,nrA,nrB,datasetType)

if strcmp(datasetType,'standard1500')
    if strcmp(section,'all')
        stem = 'fluxSumStatst_all';
    elseif strcmp(section,'allrank')
        stem = 'fluxSumStatsallrank';
    else
        stem = ['fluxSumStats_' section];
    end

    filename = [modelPrefix '_' stem '_' treatment '_4_' num2str(nrA) ...
        '_VS_' num2str(nrB) '1500.xlsx'];

elseif strcmp(datasetType,'onoff1000')
    if strcmp(section,'all')
        stem = 'fluxSumStatst_all';
    elseif strcmp(section,'allrank')
        stem = 'fluxSumStatsallrank';
    elseif strcmp(treatment,'medium')
        if strcmp(section,'Glycolysis')
            stem = 'fluxSumStatst_Glyco';
        elseif strcmp(section,'PPP')
            stem = 'fluxSumStatst_PPP';
        elseif strcmp(section,'Pyruvate')
            stem = 'fluxSumStatst_Pyr';
        elseif strcmp(section,'Glutathione')
            stem = 'fluxSumStatst_Glutath';
        else
            stem = ['fluxSumStats_' section];
        end
    else
        stem = ['fluxSumStats_' section];
    end

    filename = [modelPrefix '_' stem '_' treatment '_4_' num2str(nrA) ...
        '_vs_' num2str(nrB) '.xlsx'];
end

file = fullfile(outputFolder,filename);

end


function deleteIfPresentLegacy(file)
if exist(file,'file')
    delete(file);
end
end
