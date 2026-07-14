%% dotplots_combined_paper.m
% MATLAB R2019b-compatible.
% Creates the three combined paper figures:
%   d) Glycolysis metabolites
%   e) Pyruvate-metabolism metabolites
%   f) Pyruvate-metabolism reactions
%
% Reads stored fc/log2FC values directly from Excel outputs.
% Does not recalculate fold change, SNR, or p-values.

clear;
clc;
close all;

%% Relative folders
fluxSumStandardFolder = '../data/flux_sum/standard';
fluxSumOnOffFolder = '../data/flux_sum/onoff';
fluxSamplingStandardFolder = '../data/flux_sampling/standard';
fluxSamplingOnOffFolder = '../data/flux_sampling/onoff';
outputFolder = '../results/plots/combined';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%% Input files and condition order
fluxSumFiles = {
    fullfile(fluxSumOnOffFolder,'LIHC_fluxSumStatst_all_medium_4_7_VS_8.xlsx')
    fullfile(fluxSumOnOffFolder,'KIRC_fluxSumStatst_all_medium_4_1_VS_2.xlsx')
    fullfile(fluxSumStandardFolder,'Huh7_fluxSumStatst_all_medium_4_5_VS_71500.xlsx')
    fullfile(fluxSumStandardFolder,'769-P_fluxSumStatst_all_medium_4_1_VS_31500.xlsx')
};

samplingFiles = {
    fullfile(fluxSamplingOnOffFolder,'LIHC-ONOFF_FluxSamplingStatst_all_medium_5onoff7_vs_8.xlsx')
    fullfile(fluxSamplingOnOffFolder,'KIRC-ONOFF_FluxSamplingStatst_all_medium_5onoff1_vs_2.xlsx')
    fullfile(fluxSamplingStandardFolder,'Huh7_FluxSamplingStatst_all_medium_55_vs_7_1500.xlsx')
    fullfile(fluxSamplingStandardFolder,'769-P_FluxSamplingStatst_all_medium_51_vs_3_1500.xlsx')
};

conditions = {'LIHC';'KIRC';'Huh7';'769-P'};

%% Figure d: Glycolysis metabolites
glycolysisIDs = {
    'glc_D[c]';'g6p[c]';'f6p[c]';'fdp[c]';'g3p[c]';'dhap[c]'; ...
    '13dpg[c]';'3pg[c]';'2pg[c]';'pep[c]';'pyr[c]';'lac_L[e]'
};

glycolysisLabels = {
    'glc_D';'g6p';'f6p';'fdp';'g3p';'dhap'; ...
    '1,3dpg';'3pg';'2pg';'pep';'pyr';'lac_L'
};

glycolysisNames = {
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

glycolysisFC = readSelectedFeatures(fluxSumFiles,glycolysisIDs,'',false);

makePaperDotPlot(glycolysisFC,glycolysisLabels,conditions, ...
    'Glycolysis Pathway', ...
    fullfile(outputFolder,'Figure_d_Glycolysis_Pathway'),true);

writeSelectedMatrix(fullfile(outputFolder,'Figure_d_Glycolysis_data.xlsx'), ...
    glycolysisIDs,glycolysisNames,conditions,glycolysisFC);

%% Figure e: Pyruvate-metabolism metabolites
pyruvateMetIDs = {
    'mthgxl[c]';'12ppd_S[c]';'12ppd_R[c]';'lald_L[c]'; ...
    'lald_L[m]';'lald_D[c]';'ac[m]';'aac[c]';'acetol[c]'; ...
    'lgt_S[c]';'lgt_S[m]';'lac_D[c]';'lac_D[m]';'lac_L[m]'; ...
    'pyr[m]';'oaa[m]';'acacoa[c]';'acacoa[m]';'gthrd[c]';'gthrd[m]'
};

pyruvateMetLabels = {
    'mthgxl (c)';'12ppd_S (c)';'12ppd_R (c)';'lald_L (c)'; ...
    'lald_L (m)';'lald_D (c)';'ac (m)';'aac (c)';'acetol (c)'; ...
    'lgt_S (c)';'lgt_S (m)';'lac_D (c)';'lac_D (m)'; ...
    'lac_L (m)';'pyr (m)';'oaa (m)';'acacoa (c)'; ...
    'acacoa (m)';'gthrd (c)';'gthrd (m)'
};

pyruvateFC = readSelectedFeatures(fluxSumFiles,pyruvateMetIDs,'',false);

makePaperDotPlot(pyruvateFC,pyruvateMetLabels,conditions, ...
    'Pyruvate metabolism Pathway', ...
    fullfile(outputFolder,'Figure_e_Pyruvate_Metabolites'),true);

writeSelectedMatrix(fullfile(outputFolder,'Figure_e_Pyruvate_metabolites_data.xlsx'), ...
    pyruvateMetIDs,pyruvateMetLabels,conditions,pyruvateFC);

%% Figure f: Pyruvate-metabolism reactions
selectedReactions = {
    'MGSA';'MGSA2';'ALR2';'ALR3';'ALCD21_L';'ALCD22_L'; ...
    'LCADI';'GLYOX';'LCADI_D';'HMR_8501';'ALCD22_D'; ...
    'ALCD21_D';'LCADIm';'LGTHL';'PPDOy';'ME2';'LDH_Lm';'HMR_3855'
};

reactionFC = readSelectedFeatures( ...
    samplingFiles,selectedReactions,'Pyruvate metabolism',true);

makePaperDotPlot(reactionFC,selectedReactions,conditions, ...
    'Pyruvate metabolism Pathway (reactions)', ...
    fullfile(outputFolder,'Figure_f_Pyruvate_Reactions'),true);

writeSelectedMatrix(fullfile(outputFolder,'Figure_f_Pyruvate_reactions_data.xlsx'), ...
    selectedReactions,selectedReactions,conditions,reactionFC);

fprintf('\nFinished combined figures.\n');
fprintf('Output folder:\n%s\n',outputFolder);

%% Local functions
function result = readSelectedFeatures(files,selectedIDs,pathwayName,useSubsystem)
result = nan(numel(files),numel(selectedIDs));
for f = 1:numel(files)
    filename = files{f};
    if ~exist(filename,'file')
        error('File not found:\n%s',filename);
    end
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

function makePaperDotPlot(values,xLabels,conditions,plotTitle,outputStub,showGroups)
nConditions = size(values,1);
nFeatures = size(values,2);
colorLimit = 2;
minMarkerArea = 28;
maxMarkerArea = 420;
sizeExponent = 1.10;
if nFeatures <= 12
    pixelsPerFeature = 72; labelFontSize = 12; labelRotation = 50;
elseif nFeatures <= 18
    pixelsPerFeature = 68; labelFontSize = 11; labelRotation = 55;
else
    pixelsPerFeature = 72; labelFontSize = 10; labelRotation = 60;
end
figWidth = max(1100,pixelsPerFeature*nFeatures+420);
fig = figure('Position',[60 60 figWidth 780],'Color','w','PaperPositionMode','auto');
ax = axes('Parent',fig,'Position',[0.145 0.31 0.635 0.56]);
hold(ax,'on');
colormap(ax,blueGrayRed(256));
caxis(ax,[-colorLimit colorLimit]);
if nConditions == 4
    yPositions = [4.2 3.2 1.2 0.2];
    yTicks = [0.2 1.2 3.2 4.2];
    yLabels = conditions(end:-1:1);
    yLimits = [-0.35 4.75];
else
    yPositions = nConditions:-1:1;
    yTicks = 1:nConditions;
    yLabels = conditions(end:-1:1);
    yLimits = [0.5 nConditions+0.5];
end
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
set(ax,'XLim',[0.35 nFeatures+0.65],'YLim',yLimits, ...
    'XTick',1:nFeatures,'XTickLabel',xLabels,'XTickLabelRotation',labelRotation, ...
    'YTick',yTicks,'YTickLabel',yLabels,'TickLabelInterpreter','none', ...
    'FontName','Arial','FontSize',labelFontSize,'FontWeight','bold', ...
    'TickDir','out','Box','off','LineWidth',0.8);
title(ax,plotTitle,'FontName','Arial','FontSize',18,'FontWeight','normal');
ax.XColor = [0.15 0.15 0.15];
ax.YColor = [0.15 0.15 0.15];
if showGroups && nConditions == 4
    bracketX = -1.25; textX = -1.72;
    line(ax,[bracketX bracketX],[3.0 4.45],'Color',[0.15 0.15 0.15],'LineWidth',1.3,'Clipping','off');
    line(ax,[bracketX bracketX],[0.0 1.45],'Color',[0.15 0.15 0.15],'LineWidth',1.3,'Clipping','off');
    text(ax,textX,3.72,'Patients','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Arial','FontSize',15,'Clipping','off');
    text(ax,textX,0.72,'Cell lines','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Arial','FontSize',15,'Clipping','off');
end
cb = colorbar(ax);
set(cb,'Position',[0.805 0.59 0.024 0.245],'Ticks',[-2 -1 0 1 2], ...
    'TickLabels',{'<-2','-1','0','1','>2'},'FontName','Arial','FontSize',10);
title(cb,'log2\_FC','Interpreter','tex','FontName','Arial','FontSize',11,'FontWeight','normal');
print(fig,[outputStub '.svg'],'-dsvg');
print(fig,[outputStub '.pdf'],'-dpdf','-bestfit');
close(fig);
end

function [headers,raw] = readExcelRaw(filename)
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
