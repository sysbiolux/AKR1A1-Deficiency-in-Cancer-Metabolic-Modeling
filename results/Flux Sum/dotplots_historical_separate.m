%% dotplots_historical_separate.m
% MATLAB R2019b-compatible.
%
% Creates separate historical figure families represented by the SVG files.
% Reads stored fc/log2FC values directly from the analysis workbooks.
% Does not recalculate FC, SNR, p-values, or pathway statistics.

clear;
clc;
close all;

%% Relative folders
fluxSumStandardFolder = '../data/flux_sum/standard';
fluxSamplingStandardFolder = '../data/flux_sampling/standard';
fluxSamplingOnOffFolder = '../data/flux_sampling/onoff';
outputStandardFolder = '../results/plots/historical/standard1500';
outputOnOffFolder = '../results/plots/historical/onoff1000';

if ~exist(outputStandardFolder,'dir'), mkdir(outputStandardFolder); end
if ~exist(outputOnOffFolder,'dir'), mkdir(outputOnOffFolder); end

%% Standard 1500 comparison files
standardFluxSumFiles = {
    fullfile(fluxSumStandardFolder,'769-P_fluxSumStatst_all_medium_4_1_VS_21500.xlsx')
    fullfile(fluxSumStandardFolder,'769-P_fluxSumStatst_all_medium_4_1_VS_41500.xlsx')
    fullfile(fluxSumStandardFolder,'Huh7_fluxSumStatst_all_medium_4_5_VS_61500.xlsx')
    fullfile(fluxSumStandardFolder,'Huh7_fluxSumStatst_all_medium_4_5_VS_81500.xlsx')
};

standardFluxSamplingFiles = {
    fullfile(fluxSamplingStandardFolder,'769-P_FluxSamplingStatst_all_medium_51_vs_2_1500.xlsx')
    fullfile(fluxSamplingStandardFolder,'769-P_FluxSamplingStatst_all_medium_51_vs_4_1500.xlsx')
    fullfile(fluxSamplingStandardFolder,'Huh7_FluxSamplingStatst_all_medium_55_vs_6_1500.xlsx')
    fullfile(fluxSamplingStandardFolder,'Huh7_FluxSamplingStatst_all_medium_55_vs_8_1500.xlsx')
};

standardConditions = {'RCCsi1';'RCCsi2';'HCCsi1';'HCCsi2'};

%% TCGA ON/OFF comparison files
onoffFluxSamplingFiles = {
    fullfile(fluxSamplingOnOffFolder,'KIRC-ONOFF_FluxSamplingStatst_all_medium_5onoff1_vs_2.xlsx')
    fullfile(fluxSamplingOnOffFolder,'KIRP-ONOFF_FluxSamplingStatst_all_medium_5onoff3_vs_4.xlsx')
    fullfile(fluxSamplingOnOffFolder,'KICH-ONOFF_FluxSamplingStatst_all_medium_5onoff5_vs_6.xlsx')
    fullfile(fluxSamplingOnOffFolder,'LIHC-ONOFF_FluxSamplingStatst_all_medium_5onoff7_vs_8.xlsx')
};

onoffConditions = {'KIRC';'KIRP';'KICH';'LIHC'};

%% Explicit metabolite order for historical standard plots
glycolysisMetIDs = {
    'glc_D[c]';'g6p[c]';'f6p[c]';'fdp[c]';'g3p[c]';'dhap[c]'; ...
    '13dpg[c]';'3pg[c]';'2pg[c]';'pep[c]';'pyr[c]';'lac_L[e]'
};

glycolysisMetNames = {
    'D-Glucose'
    'Glucose 6-phosphate'
    'Fructose 6-phosphate'
    'Fructose 1,6-bisphosphate'
    'Glyceraldehyde 3-phosphate'
    'Dihydroxyacetone phosphate'
    '1,3-Bisphosphoglycerate'
    '3-Phosphoglycerate'
    '2-Phosphoglycerate'
    'Phosphoenolpyruvate'
    'Pyruvate'
    'L-Lactate'
};

pyruvateMetIDs = {
    'mthgxl[c]';'12ppd_S[c]';'12ppd_R[c]';'lald_L[c]'; ...
    'lald_L[m]';'lald_D[c]';'ac[m]';'aac[c]';'acetol[c]'; ...
    'lgt_S[c]';'lgt_S[m]';'lac_D[c]';'lac_D[m]';'lac_L[m]'; ...
    'pyr[m]';'oaa[m]';'acacoa[c]';'acacoa[m]';'gthrd[c]';'gthrd[m]'
};

pyruvateMetNames = pyruvateMetIDs;

%% Standard 1500 metabolite figures
glycolysisMetFC = readSelectedFeatures(standardFluxSumFiles,glycolysisMetIDs,'',false);
makeHistoricalDotPlot(glycolysisMetFC,glycolysisMetIDs,standardConditions, ...
    'Fluxsum: Glycolysis-gluconeogenesis pathway', ...
    fullfile(outputStandardFolder,'Flux_Sampling_Glycolysis-Gluconeogenesis_Pathway1500_MetsIDsby1to2'));
makeHistoricalDotPlot(glycolysisMetFC,glycolysisMetNames,standardConditions, ...
    'Fluxsum: Glycolysis-gluconeogenesis pathway', ...
    fullfile(outputStandardFolder,'Flux_Sampling_Glycolysis-Gluconeogenesis_Pathway1500_MetsNamesby1to2'));

pyruvateMetFC = readSelectedFeatures(standardFluxSumFiles,pyruvateMetIDs,'',false);
makeHistoricalDotPlot(pyruvateMetFC,pyruvateMetIDs,standardConditions, ...
    'Fluxsum: Pyruvate metabolism pathway', ...
    fullfile(outputStandardFolder,'Flux_Sampling_Pyruvate_Metabolism_Pathway1500_MetsIDsby1to2'));
makeHistoricalDotPlot(pyruvateMetFC,pyruvateMetNames,standardConditions, ...
    'Fluxsum: Pyruvate metabolism pathway', ...
    fullfile(outputStandardFolder,'Flux_Sampling_Pyruvate_Metabolism_Pathway1500_MetsNamesby1to2'));

%% Standard 1500 reaction figures
makeReactionFigurePair(standardFluxSamplingFiles,standardConditions, ...
    'Glycolysis/gluconeogenesis','Glycolysis-Gluconeogenesis', ...
    '1500',outputStandardFolder,true);
makeReactionFigurePair(standardFluxSamplingFiles,standardConditions, ...
    'Pyruvate metabolism','Pyruvate_Metabolism', ...
    '1500',outputStandardFolder,true);

%% TCGA ON/OFF reaction figures
makeReactionFigurePair(onoffFluxSamplingFiles,onoffConditions, ...
    'Glycolysis/gluconeogenesis','Glycolysis-Gluconeogenesis', ...
    'TCGA',outputOnOffFolder,false);
makeReactionFigurePair(onoffFluxSamplingFiles,onoffConditions, ...
    'Pyruvate metabolism','Pyruvate_Metabolism', ...
    'TCGA',outputOnOffFolder,false);

fprintf('\nFinished historical separate plots.\n');

%% Local functions
function makeReactionFigurePair(files,conditions,pathwayName,pathwayFileLabel,datasetLabel,outputFolder,standardUnderscore)
[reactionIDs,reactionNames,values] = readAllPathwayReactions(files,pathwayName);
if standardUnderscore
    idsStub = ['Flux_Sampling_' pathwayFileLabel '_Pathway' datasetLabel '_ReactionsIDsby1to2_'];
    namesStub = ['Flux_Sampling_' pathwayFileLabel '_Pathway' datasetLabel '_ReactionsNamesby1to2_'];
else
    idsStub = ['Flux_Sampling_' pathwayFileLabel '_Pathway' datasetLabel '_ReactionsIDsby1to2'];
    namesStub = ['Flux_Sampling_' pathwayFileLabel '_Pathway' datasetLabel '_ReactionsNamesby1to2'];
end
makeHistoricalDotPlot(values,reactionIDs,conditions,['Flux sampling: ' pathwayName],fullfile(outputFolder,idsStub));
makeHistoricalDotPlot(values,reactionNames,conditions,['Flux sampling: ' pathwayName],fullfile(outputFolder,namesStub));
writeSelectedMatrix(fullfile(outputFolder,[idsStub '_data.xlsx']),reactionIDs,reactionNames,conditions,values);
end

function [featureIDs,featureNames,values] = readAllPathwayReactions(files,pathwayName)
allIDs = {};
allNames = {};
for f = 1:numel(files)
    [headers,raw] = readExcelRaw(files{f});
    idCol = findColumn(headers,{'Row','rxns','rxnID','ReactionID','...1'});
    nameCol = findColumn(headers,{'rxnNames','ReactionName','reactionNames'});
    subsystemCol = findColumn(headers,{'subSystems','subSystem','subsys','Subsystem'});
    if isempty(idCol), idCol = 1; end
    if isempty(subsystemCol), error('No subsystem column found in:\n%s',files{f}); end
    for r = 2:size(raw,1)
        currentSubsystem = strtrim(cellToText(raw{r,subsystemCol}));
        if ~strcmpi(currentSubsystem,pathwayName), continue; end
        currentID = strtrim(cellToText(raw{r,idCol}));
        if isempty(currentID), continue; end
        if isempty(find(strcmpi(allIDs,currentID),1))
            allIDs{end+1,1} = currentID; %#ok<AGROW>
            if isempty(nameCol)
                allNames{end+1,1} = currentID; %#ok<AGROW>
            else
                currentName = strtrim(cellToText(raw{r,nameCol}));
                if isempty(currentName), currentName = currentID; end
                allNames{end+1,1} = currentName; %#ok<AGROW>
            end
        end
    end
end
featureIDs = allIDs;
featureNames = allNames;
values = readSelectedFeatures(files,featureIDs,pathwayName,true);
end

function result = readSelectedFeatures(files,selectedIDs,pathwayName,useSubsystem)
result = nan(numel(files),numel(selectedIDs));
for f = 1:numel(files)
    filename = files{f};
    if ~exist(filename,'file'), error('File not found:\n%s',filename); end
    [headers,raw] = readExcelRaw(filename);
    idCol = findColumn(headers,{'Row','mets','Metabolite_ID','MetaboliteID','metID','rxns','rxnID','ReactionID','...1'});
    fcCol = findColumn(headers,{'log2FC','log2_FC','log2.fc','fc'});
    if isempty(idCol), idCol = 1; end
    if isempty(fcCol), error('No fc/log2FC column found in:\n%s',filename); end
    if useSubsystem
        subsystemCol = findColumn(headers,{'subSystems','subSystem','subsys','Subsystem'});
        if isempty(subsystemCol), error('No subsystem column found in:\n%s',filename); end
    else
        subsystemCol = [];
    end
    for j = 1:numel(selectedIDs)
        rowIndex = [];
        for r = 2:size(raw,1)
            currentID = strtrim(cellToText(raw{r,idCol}));
            idMatches = strcmpi(currentID,strtrim(selectedIDs{j})) || ...
                strcmpi(normalizeMetaboliteID(currentID),normalizeMetaboliteID(selectedIDs{j}));
            if useSubsystem
                currentSubsystem = strtrim(cellToText(raw{r,subsystemCol}));
                pathwayMatches = strcmpi(currentSubsystem,pathwayName);
            else
                pathwayMatches = true;
            end
            if idMatches && pathwayMatches
                rowIndex = r;
                break;
            end
        end
        if isempty(rowIndex)
            warning('Feature not found in %s: %s',filename,selectedIDs{j});
        else
            result(f,j) = cellToNumber(raw{rowIndex,fcCol});
        end
    end
end
end

function makeHistoricalDotPlot(values,xLabels,conditions,plotTitle,outputStub)
nConditions = size(values,1);
nFeatures = size(values,2);
colorLimit = 2;
minMarkerArea = 28;
maxMarkerArea = 420;
sizeExponent = 1.10;
if nFeatures <= 12
    pixelsPerFeature = 72; labelFontSize = 11; labelRotation = 50;
elseif nFeatures <= 20
    pixelsPerFeature = 68; labelFontSize = 10; labelRotation = 55;
else
    pixelsPerFeature = 62; labelFontSize = 9; labelRotation = 60;
end
figWidth = max(1100,pixelsPerFeature*nFeatures+350);
fig = figure('Position',[60 60 figWidth 700],'Color','w','PaperPositionMode','auto');
ax = axes('Parent',fig,'Position',[0.11 0.30 0.72 0.58]);
hold(ax,'on');
colormap(ax,blueGrayRed(256));
caxis(ax,[-colorLimit colorLimit]);
yPositions = nConditions:-1:1;
for row = 1:nConditions
    for col = 1:nFeatures
        value = values(row,col);
        if isnan(value), continue; end
        clippedValue = max(min(value,colorLimit),-colorLimit);
        magnitude = abs(clippedValue)/colorLimit;
        markerArea = minMarkerArea + (maxMarkerArea-minMarkerArea)*magnitude^sizeExponent;
        scatter(ax,col,yPositions(row),markerArea,clippedValue,'filled', ...
            'MarkerEdgeColor',[0.70 0.70 0.70],'LineWidth',0.30);
    end
end
set(ax,'XLim',[0.35 nFeatures+0.65],'YLim',[0.5 nConditions+0.5], ...
    'XTick',1:nFeatures,'XTickLabel',xLabels,'XTickLabelRotation',labelRotation, ...
    'YTick',1:nConditions,'YTickLabel',conditions(end:-1:1), ...
    'TickLabelInterpreter','none','FontName','Arial','FontSize',labelFontSize, ...
    'FontWeight','bold','TickDir','out','Box','off','LineWidth',0.8);
title(ax,plotTitle,'FontName','Arial','FontSize',16,'FontWeight','normal');
cb = colorbar(ax);
set(cb,'Ticks',[-2 -1 0 1 2],'TickLabels',{'<-2','-1','0','1','>2'},'FontName','Arial','FontSize',10);
title(cb,'log2\_FC','Interpreter','tex','FontName','Arial','FontSize',10,'FontWeight','normal');
print(fig,[outputStub '.svg'],'-dsvg');
print(fig,[outputStub '.pdf'],'-dpdf','-bestfit');
close(fig);
end

function [headers,raw] = readExcelRaw(filename)
if ~exist(filename,'file'), error('Excel file not found:\n%s',filename); end
[~,baseName,extension] = fileparts(filename);
localCopy = fullfile(tempdir,[baseName '_' datestr(now,'yyyymmdd_HHMMSSFFF') extension]);
[copyOK,copyMessage] = copyfile(filename,localCopy,'f');
if ~copyOK, error('Could not copy workbook:\n%s\n%s',filename,copyMessage); end
cleanupObject = onCleanup(@() deleteTemporaryWorkbook(localCopy));
[~,~,raw] = xlsread(localCopy,'','','basic');
clear cleanupObject;
if isempty(raw), error('Workbook is empty or unreadable:\n%s',filename); end
headers = raw(1,:);
end

function deleteTemporaryWorkbook(filename)
if exist(filename,'file')
    try, delete(filename); catch, end
end
end

function columnIndex = findColumn(headers,possibleNames)
columnIndex = [];
for c = 1:numel(headers)
    currentHeader = strtrim(cellToText(headers{c}));
    for p = 1:numel(possibleNames)
        if strcmpi(currentHeader,possibleNames{p})
            columnIndex = c;
            return;
        end
    end
end
end

function normalizedID = normalizeMetaboliteID(value)
normalizedID = lower(strtrim(cellToText(value)));
normalizedID = regexprep(normalizedID,'\s+','');
normalizedID = regexprep(normalizedID,'\(([a-z])\)$','[$1]');
normalizedID = regexprep(normalizedID,'_([a-z])$','[$1]');
end

function textValue = cellToText(value)
if isempty(value), textValue = '';
elseif ischar(value), textValue = value;
elseif isstring(value), textValue = char(value);
elseif isnumeric(value) && isscalar(value), textValue = num2str(value);
else, textValue = '';
end
end

function numberValue = cellToNumber(value)
if isnumeric(value) && isscalar(value), numberValue = double(value);
elseif ischar(value), numberValue = str2double(strrep(value,',','.'));
elseif isstring(value), numberValue = str2double(strrep(char(value),',','.'));
else, numberValue = NaN;
end
end

function cmap = blueGrayRed(n)
if nargin < 1, n = 256; end
blue = [0.20 0.31 0.65]; gray = [0.78 0.78 0.78]; red = [0.93 0.29 0.20];
n1 = floor(n/2); n2 = n-n1;
firstHalf = [linspace(blue(1),gray(1),n1)',linspace(blue(2),gray(2),n1)',linspace(blue(3),gray(3),n1)'];
secondHalf = [linspace(gray(1),red(1),n2)',linspace(gray(2),red(2),n2)',linspace(gray(3),red(3),n2)'];
cmap = [firstHalf;secondHalf];
end

function writeSelectedMatrix(filename,ids,names,conditions,values)
T = table(ids(:),names(:),'VariableNames',{'ID','Name'});
for c = 1:numel(conditions)
    validName = makeSafeVariableName(conditions{c});
    T.(validName) = values(c,:)';
end
if exist(filename,'file'), delete(filename); end
writetable(T,filename);
end

function validName = makeSafeVariableName(inputName)
validName = char(inputName);
validName = strtrim(validName);
validName = regexprep(validName,'[^A-Za-z0-9_]','_');
if isempty(validName), validName = 'Variable'; end
if ~isletter(validName(1)), validName = ['x_' validName]; end
if length(validName) > namelengthmax, validName = validName(1:namelengthmax); end
end
