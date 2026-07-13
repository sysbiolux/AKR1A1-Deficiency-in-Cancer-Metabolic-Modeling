clearvars -except solverOK, clc, close all
%%
%sc12 consensus 769-P
nrA=1;
nrB=3;
model_prefix={'769-P'};
model_prefix_str = char(model_prefix);  % Convert to string if it's a single cell
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type 7
%%
%sc1 condition 1
nrA=1;
nrB=2;
model_prefix={'769-P'};
model_prefix_str = char(model_prefix);  % Convert to string if it's a single cell
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type 7
%%
%sc2 condition 2
nrA=1;
nrB=4;
model_prefix={'769-P'};
model_prefix_str = char(model_prefix);  % Convert to string if it's a single cell
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type 7
%%
%sc12 consensus Huh7
nrA=5;
nrB=7;
model_prefix={'Huh7'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
%sc1 condition 1
nrA=5;
nrB=6;
model_prefix={'Huh7'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
%sc2 condition 2
nrA=5;
nrB=8;
model_prefix={'Huh7'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
%'model_KIRC_ON'
%'model_KIRC_OFF'
nrA=1;
nrB=2;
model_prefix={'KIRC'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type KIRC
%%
%'model_KIRP_ON'
%'model_KIRP_OFF'
nrA=3;
nrB=4;
model_prefix={'KIRP'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
%'model_KICH_ON'
%'model_KICH_OFF'
nrA=5;
nrB=6;
model_prefix={'KICH'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
%'model_LIHC_ON'
%'model_LIHC_OFF'
nrA=7;
nrB=8;
model_prefix={'LIHC'};
model_prefix_str = char(model_prefix);
analyzeFluxSumData(model_prefix_str, nrA, nrB);  % For cell type H
%%
function analyzeFluxSumData(model_prefix_str, nrA, nrB)
% Clear the workspace except for necessary variables, close figures, and clear command window
clearvars -except model_prefix_str nrA nrB
clc;
close all;

%load('../consistent_model.mat')
load('consistent_model.mat');
model_orig=consistent_model;
% model=model_orig;

epsilon=1e-4
metsOfInterest=model_orig.mets;
metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));

% Define treatments and file path template
treatments = {'medium', 'KO_GLO1'};
%basePath = '../path/to/data/'; % Adjust as necessary

% Iterate over treatments
for t = 1:length(treatments)
    treatment = treatments{t};
    disp(['Processing: ' model_prefix_str ' - ' treatment]);
  
    % Build file paths dynamically based on treatment and modelPrefix
    %file=['SamplingResults_medium_1500_model_' num2str(nrA) '.mat'];
    %file=['samplingResults_medium_1000_model_' num2str(nrA) '.mat'];
    %'model_KIRC_ON'
    %file=['samplingResults_mediumonoff_1000_model_' num2str(nrA) '.mat'];
    file=['SamplingResults_medium_1500_model_' num2str(nrA) '.mat'];
    % Load data for control model
    load(file, 'x');
    model1 = x.modelSampling;
    data1 = x.samples(:, 1:500);
    size(data1)
    
    % Load data for condition model
    %file=['SamplingResults_medium_1500_model_' num2str(nrB) '.mat'];
    %file=['samplingResults_medium_1000_model_' num2str(nrB) '.mat'];
    %'model_KIRC_ON'
    %file=['samplingResults_mediumonoff_1000_model_' num2str(nrB) '.mat'];
    file=['SamplingResults_medium_1500_model_' num2str(nrB) '.mat'];
    load(file, 'x');
    model2=x.modelSampling;
    data2=x.samples(:,1:500);
    size(data2)
    
    disp('... model & data loading done ...')
    
    %% calculate flux sum per metabolite
    res1=[];
    for counter=1:size(data1,2)
        v=data1(:,counter);
        temp=repmat(v',size(model1.S,1),1);
        fluxes=model1.S.*temp;
        fluxSumP=full(sum((fluxes>0).*fluxes,2));
        fluxSumN=full(sum((fluxes<0).*fluxes,2));
        res1=[res1, fluxSumP];
    end
    disp('... fluxSum A calculated ...')
    
    res2=[];
    for counter=1:size(data2,2)
        v=data2(:,counter);
        temp=repmat(v',size(model2.S,1),1);
        fluxes=model2.S.*temp;
        fluxSumP=full(sum((fluxes>0).*fluxes,2));
        fluxSumN=full(sum((fluxes<0).*fluxes,2));
        res2=[res2, fluxSumP];
    end
    disp('... fluxSum B calculated ...')
    
    figure
    boxplot(res1(1:30,:)','Labels',model1.mets(1:30))
    set(gca,'FontSize',10,'XTickLabelRotation',45)
    
    mean(res1(find(ismember(model1.mets,'glc_D[c]')),:))
    mean(res2(find(ismember(model2.mets,'glc_D[c]')),:))
    
    %% mapping to original model
    [C,IA,IB] = intersect(model_orig.mets,model1.mets,'stable');
    resA=zeros(numel(model_orig.mets),500);
    resA(IA,:)=res1;
    
    [C,IA,IB] = intersect(model_orig.mets,model2.mets,'stable');
    resB=zeros(numel(model_orig.mets),500);
    resB(IA,:)=res2;
    
    mean(resA(find(ismember(model_orig.mets,'glc_D[c]')),:))
    mean(resB(find(ismember(model_orig.mets,'glc_D[c]')),:))
    
    %contains the common metabolite IDs.
    %IA contains the indices of the common metabolites in model_orig
    %IB contains the indices of the common metabolites in model1.
    
    %% statistical test
    stats=[];
    
    for counter=1:size(resA,1)
        A=resA(counter,:);
        B=resB(counter,:);
        
        P=ranksum(A,B);
        
        %stats=[stats; mean(A), mean(B), mean(B)-mean(A), mean(B)/mean(A), P -log10(P)];
        % Check if the mean of B is zero
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
        snr =(mean(B1) - mean(A1)) / (std(A) + std(B))
        log2FC=log2(abs(mean(B1)/mean(A1)));
        stats=[stats; mean(A), mean(B), std(A), std(B), mean(B)-mean(A), log2FC, snr, P, -log2(P)];
    end
    stats(1:10,:)
    
    figure
    histogram(stats(:,9))
    %title('P values (-log10)')
    title([model_prefix_str 'P values (-log2)'])
    figure
    %hist(log10(stats(:,4)))
    histogram(stats(:,6))
    %title('log10 foldchange (mean(B)/mean(A))')
    title([model_prefix_str 'log2 foldchange (mean(B)/mean(A))'])
    
    figure
    % plot(log10(stats(:,4)),stats(:,6),'*')
    % title('vulcano: log10 foldchange vs -log10(P)')
    plot(stats(:,6),stats(:,9),'*')
    title([model_prefix_str 'vulcano: log2 foldchange vs -log2(P)'])
    
    %% writetable in excel
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    %dn=find((log10(stats(:,4))<0).*(stats(:,6)>0).*(stats(:,3)<-cutoffFluxDiff));
    %stats=setdiff(stats,cofactors);
    statst=array2table(stats,'RowNames',model_orig.mets,'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    statst.metNames=model_orig.metNames;
    %statst=sortrows(statst,11,'descend')
    
    file=[model_prefix_str '_fluxSumStatst_all_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(statst,file,'WriteRowNames',true)
    
    %% by subSystem glycolysis
    metsOfInterest={};
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    pathway='Glycolysis/gluconeogenesis';
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    %metListselect={'glc_D[c]','g6p[c]','fdp[c]','g3p[c]','dhap[c]','3pg[c]','2pg[c]','pep[c]','pyr[c]','lac_L[e]'};
    metListidx=find(ismember(model_orig.mets, metList));
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    % upt=[upt(:,end), upt(:,1:(end-1))];
    %upt=sortrows(upt,6,'descend')
    upt=sortrows(upt,7,'descend') %by SNR
    
    % dn=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,6)<-cutoffFluxDiff));
    % dn=setdiff(dn,cofactors);
    % dn=intersect(dn,metListidx);
    % dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    % dnt.metNames=model_orig.metNames(dn);
    % dnt=sortrows(dnt,9,'descend');
    
    file=[model_prefix_str '_fluxSumStats_Glycolysis_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    
    %% PPP
    metsOfInterest={};
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    pathway='Pentose phosphate pathway';
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    %metListselect={'ru5p_D[c]','s7p[c]','e4p[c]','nadph[c]','r5p[c]','xu5p_D[c]','f6p[c]','g3p[c]'}
    metListidx=find(ismember(model_orig.mets, metList));
    
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    % upt=[upt(:,end), upt(:,1:(end-1))];
    %upt=sortrows(upt,6,'descend')
    
    % dn=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,6)<-cutoffFluxDiff));
    % dn=setdiff(dn,cofactors);
    % dn=intersect(dn,metListidx);
    % dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    % dnt.metNames=model_orig.metNames(dn);
    % dnt=sortrows(dnt,9,'descend');
    upt=sortrows(upt,7,'descend') %by SNR
    
    file=[model_prefix_str '_fluxSumStats_PPP_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    
    %% pyruvate
    metsOfInterest={};
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    pathway='Pyruvate metabolism';
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    metListidx=find(ismember(model_orig.mets, metList));
    
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    % upt=[upt(:,end), upt(:,1:(end-1))];
    %upt=sortrows(upt,6,'descend')
    upt=sortrows(upt,7,'descend') %by SNR
    
    % dn=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,6)<-cutoffFluxDiff));
    % dn=setdiff(dn,cofactors);
    % dn=intersect(dn,metListidx);
    % dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    % dnt.metNames=model_orig.metNames(dn);
    % dnt=sortrows(dnt,9,'descend');
    file=[model_prefix_str '_fluxSumStats_Pyruvate_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    
    %% glutathione
    metsOfInterest={};
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    pathway='Glutathione metabolism'
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    metListidx=find(ismember(model_orig.mets, metList));
    
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    % upt=[upt(:,end), upt(:,1:(end-1))];
    %upt=sortrows(upt,6,'descend')
    upt=sortrows(upt,7,'descend') %by SNR
    
    % dn=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,6)<-cutoffFluxDiff));
    % dn=setdiff(dn,cofactors);
    % dn=intersect(dn,metListidx);
    % dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    % dnt.metNames=model_orig.metNames(dn);
    % dnt=sortrows(dnt,9,'descend');
    
    file=[model_prefix_str '_fluxSumStats_Glutathione_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    %% ROS detoxification
    metsOfInterest={};
    pathway='ROS detoxification';
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    metListidx=find(ismember(model_orig.mets, metList));
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    upt=sortrows(upt,7,'descend') %by SNR
    
    file=[model_prefix_str '_fluxSumStats_ROS_detoxification_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    %% TCA
    
    metsOfInterest={};
    pathway='Citric acid cycle';
    metList=findMetsFromRxns(model_orig,model_orig.rxns(ismember(model_orig.subSystems,pathway)));
    metListidx=find(ismember(model_orig.mets, metList));
    cutoffFluxDiff=0
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    up=setdiff(up,cofactors);
    up=intersect(up,metListidx);
    up2=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    up2=setdiff(up2,cofactors);
    up2=intersect(up2,metListidx);
    up=union(up, up2);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    upt=sortrows(upt,7,'descend') %by SNR
    
    file=[model_prefix_str '_fluxSumStats_Citric_acid_cycle_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true)
    
    %% up/down
    cutoffFluxDiff=100
    metsOfInterest=model_orig.mets;
    metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));
    str = fileread('metsCofactors.txt');
    cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
    cofactors=find(ismember(model_orig.metNames,cofactorNames))
    
    % Finding upregulated metabolites:
    %  It finds the indices of metabolites that meet certain criteria for upregulation:
    %  The ratio of means (stats(:,4)) is greater than zero (indicating an increase in flux).
    %  The p-value (stats(:,6)) is greater than zero (indicating statistical significance).
    %  The difference in means (stats(:,3)) exceeds a predefined threshold cutoffFluxDiff.
    %  It removes cofactors and retains only metabolites of interest (metsOfInterest).
    
    up=find((stats(:,6)>0).*(stats(:,9)>0).*(stats(:,5)>cutoffFluxDiff));
    %up=find((log2(stats(:,6))>0).*(stats(:,9)>0));
    %up=setdiff(up,cofactors);
    up=intersect(up,metsOfInterest);
    upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB', 'stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    upt.metNames=model_orig.metNames(up);
    % upt=[upt(:,end), upt(:,1:(end-1))];
    %upt=sortrows(upt,6,'descend') %by foldchange
    upt=sortrows(upt,7,'descend'); %by SNR
    
    dn=find((stats(:,6)<0).*(stats(:,9)>0).*(stats(:,5)<-cutoffFluxDiff));
    %dn=find((log2(stats(:,6))<0).*(stats(:,9)>0));
    %dn=setdiff(dn,cofactors);
    dn=intersect(dn,metsOfInterest);
    dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','stdA','stdB','diff','fc','SNR','pValue','-log2_p'});
    dnt.metNames=model_orig.metNames(dn);
    dnt=sortrows(dnt,9,'descend'); %by pvalue
    
    file=[model_prefix_str '_fluxSumStatsallrank_' treatment '_4_' num2str(nrA) '_VS_' num2str(nrB) '1500.xlsx'];
    delete(file)
    writetable(upt,file,'WriteRowNames',true,'Sheet','Up')
    writetable(dnt,file,'WriteRowNames',true,'Sheet','Down')
    
end
end

