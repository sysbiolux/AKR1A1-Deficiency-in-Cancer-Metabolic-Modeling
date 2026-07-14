function [pathways, scores] = read_FVA_xlsx(filename)
%READ_FVA_XLSX Read pathway names and FVA similarity scores from Excel.
%
% Expected format:
%   Row 1    = headers
%   Column 1 = pathway/subsystem
%   Column 2 = FVA similarity score
%
% MATLAB R2019b compatible.

if exist(filename,'file') ~= 2
    error('Input file not found: %s',filename);
end

[~,~,raw] = xlsread(filename);

if isempty(raw)
    error('File is empty or could not be read: %s',filename);
end

if size(raw,2) < 2
    error('File %s must contain at least two columns.',filename);
end

header = raw(1,:);
pathwayCol = [];

for j = 1:numel(header)

    if isempty(header{j})
        continue
    end

    headerName = strtrim(char(string(header{j})));

    if strcmpi(headerName,'Pathways') || ...
       strcmpi(headerName,'Pathway') || ...
       strcmpi(headerName,'subSys') || ...
       strcmpi(headerName,'Subsystem')

        pathwayCol = j;
        break
    end
end

if isempty(pathwayCol)
    pathwayCol = 1;
end

scoreCol = 2;

pathwayCells = raw(2:end,pathwayCol);
scoreCells = raw(2:end,scoreCol);

pathways = strings(numel(pathwayCells),1);
scores = nan(numel(scoreCells),1);

for r = 1:numel(pathwayCells)

    if isempty(pathwayCells{r})
        pathways(r) = "";
    else
        pathways(r) = strtrim(string(pathwayCells{r}));
    end

    value = scoreCells{r};

    if isnumeric(value) && isscalar(value)
        scores(r) = double(value);
    elseif ischar(value) || isstring(value)
        scores(r) = str2double(string(value));
    else
        scores(r) = NaN;
    end
end

valid = ...
    isfinite(scores) & ...
    strlength(strtrim(pathways)) > 0;

pathways = pathways(valid);
scores = scores(valid);

if isempty(scores)
    error('No valid scores found in %s.',filename);
end

end