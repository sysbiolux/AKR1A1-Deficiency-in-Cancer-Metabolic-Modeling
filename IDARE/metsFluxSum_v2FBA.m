clearvars -except solverOK, clc, close all
% % % load consistent_model.mat
% % % model_orig=consistent_model;
% % % model=model_orig;
% % % 
% % % epsilon=1e-4
% % % 
% % % temp=find(model.c)
% % % model.rxns(temp)
% % % 
% % % load('medium_example.mat')
% % % medium_exchanges = cellstr(intersect(medium_example, model.mets))
% % % % find metabolites without formulas (and add C there)
% % % model.metFormulas(cellfun('isempty',model.metFormulas))=cellstr('C');
% % % model.metFormulas(strcmp(model.metFormulas,'X'))=cellstr('C');
% % % % biomass of medium comstrained model
% % % exRxns=model.rxns(findExcRxns(model));
% % % % medium_exchanges=strrep(medium_exchanges,'[e]','[csf]');
% % % temp=findRxnsFromMets(model,medium_exchanges);
% % % % convert
% % % exToKeep=intersect(exRxns,temp); %medium exchange reactions
% % % disp(numel(exToKeep))
% % % % remove uptake of carbon sources which are not in the medium from the model
% % % [model_cons] = constrain_model_rFASTCORMICS(model,medium_exchanges , [], 'biomass_maintenance', 'biomass_reaction');
% % % exToRemove=model.rxns(model.lb~=model_cons.lb | model.ub~= model_cons.ub);
% % % solution_cons = optimizeCbModel(model_cons,'max');
% % % disp(strcat(' biomass production after medium constraints: ',num2str(solution_cons.f)));
% % % % determine missing metabolites uptakes
% % % model=model_orig;
% % % model = changeObjective(model,'biomass_maintenance');
% % % C=find(ismember(model.rxns,'biomass_maintenance'));
% % % allowedInputRxns=find(ismember(model.rxns,exToKeep)); %indices of medium exchange reactions
% % % disp(' ')
% % % disp(strcat('allowed input rxns: ', num2str(numel(allowedInputRxns))));
% % % model.S(:,C)=model.S(:,C)*1000; %TS: multiple biomass coefficients to overcome numerical issues
% % % 
% % % % run fastcore for biomass only, with medium exchange reactions unpenalized
% % % A=fastcore_4_rfastcormics(C, model, epsilon, allowedInputRxns);
% % % model.rxns(A);
% % % disp(' ')
% % % disp('Reactions in minimal mode to produce biomass:')
% % % numel(A)
% % % % check if Biomass is in
% % % if ~isempty(setdiff(C,A))
% % %     error('No biomass produced')
% % % end
% % % disp('Missing medium metabolites (substrates & products):')
% % % neededUptakes=intersect(exToRemove, model.rxns(A))
% % % 
% % % % check biomass of minimal model with A only (including additional metabolites)
% % % model=model_orig;
% % % 
% % % model = changeObjective(model,'biomass_maintenance');
% % % model_min = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A))); %keep only A
% % % % Acc =  fastcc_4_rfastcormics(model_min,epsilon,0);
% % % % numel(Acc)
% % % [~,~,IB]=intersect(model_min.rxns,exRxns);
% % % exRxn_minimal=exRxns(IB);
% % % solution_min = optimizeCbModel(model_min,'max');
% % % disp(strcat('minimal solution: ', num2str(solution_min.f)))
% % % 
% % % [~, IA, IB]=intersect(model_min.rxns,exRxn_minimal); % all exchanges
% % % T=table(exRxn_minimal(IB),solution_min.x(IA));
% % % 
% % % % defined as export in S matrix
% % % [~,exports]=find(model_orig.S(:,ismember(model_orig.rxns, neededUptakes))<0);
% % % exports=neededUptakes(exports);
% % % 
% % % % defined as import in S matrix
% % % [~,imports]=find(model_orig.S(:,ismember(model_orig.rxns, neededUptakes))>0);
% % % imports=neededUptakes(imports);
% % % 
% % % % keep only importing new medium components (depending in export/import
% % % % definition)
% % % [~, IA, IB]=intersect(model_min.rxns, exports);
% % % exports=exports(IB(solution_min.x(IA)<0));
% % % [I, IA, IB]=intersect(model_min.rxns, exports);
% % % imports=imports(IB(solution_min.x(IA)>0));
% % % 
% % % disp('Missing medium metabolites (substrates only):')
% % % neededUptakes=union(imports, exports)
% % % 
% % % % check biomass of minimal model with medium and additional metabolites
% % % model=model_orig;
% % % model = changeObjective(model,'biomass_maintenance');
% % % [model_cons] = constrain_model_rFASTCORMICS(model,medium_exchanges , [], 'biomass_maintenance', 'biomass_maintenance');
% % % model_cons.lb(ismember(model.rxns,exToKeep))=model_orig.lb(ismember(model.rxns,exToKeep));
% % % model_cons.ub(ismember(model.rxns,exToKeep))=model_orig.ub(ismember(model.rxns,exToKeep));
% % % 
% % % model_cons.lb(ismember(model.rxns,neededUptakes))=model_orig.lb(ismember(model.rxns,neededUptakes));
% % % model_cons.ub(ismember(model.rxns,neededUptakes))=model_orig.ub(ismember(model.rxns,neededUptakes));
% % % 
% % % solution_corr = optimizeCbModel(model_cons,'max');
% % % solution_corr.f
% % % 
% % % temp=find(ismember(model_cons.rxns,'EX_glc_D[e]'));
% % % model_cons.lb(temp)=model_cons.lb(temp)*100;
% % % model_cons.ub(temp)=model_cons.ub(temp)*100;
% % % 
% % % sol=optimizeCbModel(model_cons,'max','zero')
% % % sum(abs(sol.v)>eps)
% % % v=sol.v;
% % % 
% % % temp=find(findExcRxns(model_cons));
% % % temp2=find(v(temp)<0);
% % % model_cons.rxns(temp(temp2))
% % % v(temp(temp2))

load('C:\Users\thomas.sauter\Dropbox\IDARE2\toEvelyn\modelNames_H.mat')
model=modelNames_H{1,2}
sol=optimizeCbModel(model,'max','zero')
sol.f

%% inputs
% model=model_cons
% v=sol.v %table2array(res(:,3)); %data to visualize
model=model
v=sol.x; %table2array(res(:,3)); %data to visualize


cutoffFluxSum=eps %keep metabolites above this cutoff
cutoffRxn=eps %keep rxn above this cutoff

keepSubSystemsOnly=[];
% keepSubSystemsOnly={'Glycolysis/gluconeogenesis'}
% keepSubSystemsOnly={'Glycolysis/gluconeogenesis','Pentose phosphate pathway','Citric acid cycle'}

removeCofactorsFromFile=[];
% removeCofactorsFromFile='metsCofactors.txt'

keepAllRxnsFromMets=0 %add all rxns for included metabolites, even if not in subsystem

maxNodeSize=150 %scale nodeSize according to this value
log2NodeSize=1 %1: log2(v+2), 0: no log2 scaling

outputFile='metsFluxSum_log2'

%% calculate flux sum per metabolite
temp=repmat(v',size(model.S,1),1);
fluxes=model.S.*temp;
fluxSumP=full(sum((fluxes>0).*fluxes,2));
fluxSumN=full(sum((fluxes<0).*fluxes,2));
temp=[fluxSumP, fluxSumN];

% number of overall and active reactions per metabolite
groupCount=full(sum(model.S~=0,2));
groupCountFlux=full(sum(fluxes~=0,2));

% rename metabolites for IDARE: coa[c] to M_coa__91__c__93__
temp=strcat('M_', model.mets);
temp=strrep(temp,'[','__91__');
temp=strrep(temp,']','__93__');

% nodeSize and shape and organize in table
if log2NodeSize
    nodeSize=log2(fluxSumP+2);
else
    nodeSize=fluxSumP;
end
keepMaxValue=max(nodeSize);
nodeSize=nodeSize*maxNodeSize/keepMaxValue; %scale to maxNodeSize
G=table(model.mets,temp,groupCount,groupCountFlux,fluxSumP,fluxSumN,nodeSize);
G.shape=repmat(cellstr('Elipse'),size(G,1),1);
G(1:10,:)

% top hits (largest flusSum)
Gup = sortrows(G,4,'descend');
Gdn = sortrows(G,4,'ascend');
Gup(1:30,:)
% Gdn(1:10,:)

% add keep/remove flag (keep above fluxSum cutoff)
keepremove={};
for counter=1:numel(fluxSumP)
    if fluxSumP(counter)>=cutoffFluxSum
        keepremove=[keepremove; 'keep'];
    else
        keepremove=[keepremove; 'remove'];
    end
end
keepremove(1:10)
G.keepremove=keepremove;

writetable(G(:,[2,7]),'metsFluxSum_log2',"QuoteStrings",0,'WriteVariableNames',0)

%% reactions
name=strcat('R_', model.rxns);
% R_EX_pcholn203_hs__91__e__93__
name=strrep(name,'[','__91__');
name=strrep(name,']','__93__');
name=strrep(name,'-','__45__');

% % % dico(ismember(dico(:,2),'PCHOLSTE_HSABCt'),3)= cellstr('R_PCHOLSTE_HSABCt');
% % % dico(ismember(dico(:,2),'PCHOLSTE_HSt1e'),3)= cellstr('R_PCHOLSTE_HSt1e');

if log2NodeSize
    nodeSize=log2(abs(v)+2);
else
    nodeSize=abs(v);
end
nodeSize=nodeSize*maxNodeSize/keepMaxValue; %scale to maxNodeSize (of mets)
shape={};
for counter=1:numel(v)
    if v(counter)>=cutoffRxn
        shape=[shape; 'Triangle'];
    elseif v(counter)<=-cutoffRxn
        shape=[shape; 'V'];
    else
        shape=[shape; 'none'];
    end
end
shape(1:10)

keepremove={};
for counter=1:numel(v)
    if abs(v(counter))>=cutoffRxn
        keepremove=[keepremove; 'keep'];
    else
        keepremove=[keepremove; 'remove'];
    end
end
keepremove(1:10)

out1=G(:,[2,7,8,9]);
out2=table(name, nodeSize, shape, keepremove,'VariableNames',out1.Properties.VariableNames);
out1(1:10,:)
out2(1:10,:)
out=[out1; out2];

%% subSystems only / remove cofactors / keep all rxns from mets
toKeep=[];
if ~isempty(keepSubSystemsOnly)
    for counter=1:numel(keepSubSystemsOnly)
        subsystem=keepSubSystemsOnly(counter)
        temp=find(ismember(model.subSystems,subsystem));
        temp2=findMetsFromRxns(model,model.rxns(temp)); %mets
        % rxns above cutoff only
        vs=[];
        for counter2=1:numel(temp)
            %             temp2=find(ismember(model.rxns,temp(counter2)));
            vs=[vs; v(temp(counter2))];
        end
        %        vs
        temp(abs(vs)<cutoffRxn)=[];
        
        fs=[];
        for counter2=1:numel(temp2)
            temp3=find(ismember(model.mets,temp2(counter2)));
            fs=[fs; fluxSumP(temp3)];
        end
        % fs
        temp2(fs<cutoffFluxSum)=[];
        
        temp=strcat('R_',model.rxns(temp)); %rxns to keep
        temp=strrep(temp,'[','__91__');
        temp=strrep(temp,']','__93__');
        temp=strrep(temp,'-','__45__');
        
        temp2=strcat('M_', temp2); %mets to keep
        temp2=strrep(temp2,'[','__91__');
        temp2=strrep(temp2,']','__93__');
        
        toKeep=[toKeep; unique([temp;temp2])];
    end
    temp3=find(ismember(table2cell(out(:,1)),unique(toKeep)));
    out(temp3,4)=repmat({'keep'},numel(temp3),1);
    
    temp4=setdiff(1:size(out,1),temp3);
    out(temp4,4)=repmat({'remove'},numel(temp4),1);
end

if ~isempty(removeCofactorsFromFile) %remove cofactor metabolites
    fileID = fopen(removeCofactorsFromFile);
    C = textscan(fileID,'%s','Delimiter','\n')
    C=C{:}
    fclose(fileID);
    
    temp=find(ismember(model.metNames,C));
    temp=model.mets(temp);
    temp=strcat('M_', temp);
    temp=strrep(temp,'[','__91__');
    temp=strrep(temp,']','__93__');
    temp2=find(ismember(table2cell(out(:,1)),temp));
    out(temp2,4)=repmat({'remove'},numel(temp2),1);
end

if keepAllRxnsFromMets
    temp=find(ismember(table2cell(out(:,4)),'keep'));
    temp2=find(contains(table2cell(out(:,1)),'M_'));
    temp2(temp2>numel(model.mets))=[]; %no rxns here
    temp3=intersect(temp,temp2);
    temp4=cellstr(out{temp3,1});
    temp5=eraseBetween(temp4,1,2);
    temp5=strrep(temp5,'__91__','[');
    temp5=strrep(temp5,'__93__',']') %all mets
    
    temp=findRxnsFromMets(model,temp5);  %rxns to keep
    %above cutoff only
    vs=[];
    for counter=1:numel(temp)
        temp2=find(ismember(model.rxns,temp(counter)));
        vs=[vs; v(temp2)];
    end
    temp(abs(vs)<cutoffRxn)=[];
    
    temp=strcat('R_',temp);
    temp=strrep(temp,'[','__91__');
    temp=strrep(temp,']','__93__');
    temp=strrep(temp,'-','__45__');
    
    %     temp2=strcat('M_', temp2); %mets to keep
    %     temp2=strrep(temp2,'[','__91__');
    %     temp2=strrep(temp2,']','__93__');
    
    temp3=find(ismember(table2cell(out(:,1)),temp));
    out(temp3,4)=repmat({'keep'},numel(temp3),1);
end

% final stats
temp=find(ismember(table2cell(out(:,4)),'keep'));
temp2=find(contains(table2cell(out(:,1)),'M_'));
temp2(temp2>numel(model.mets))=[]; %no rxns here
temp3=intersect(temp,temp2);
disp('')
disp('Kept metabolites:')
numel(temp3)
temp2=find(contains(table2cell(out(:,1)),'R_'));
temp3=intersect(temp,temp2);
disp('')
disp('Kept reactions:')
numel(temp3)

writetable(out,outputFile,"QuoteStrings",0,'WriteVariableNames',0)

%%
temp=find(ismember(model.mets,'val_L[e]'))
fluxes(temp,:)
temp2=find(fluxes(temp,:));
model.rxns(temp2)
printRxnFormula(model,model.rxns(temp2))

temp3=find(contains(model.rxns, 'EX_'))
v(temp3)
model.rxns(temp3)
% printRxnFormula(model,model.rxns(temp3))
table(model.rxns(temp3),v(temp3))
