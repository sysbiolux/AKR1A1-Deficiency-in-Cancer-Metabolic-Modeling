function [Table_ON, Table_OFF, colnamesON, colnamesOFF] = ...
    high_low(data, genes, colnames, geneName)

% Convert genes to a column cell array
if ischar(genes)
    genes = cellstr(genes);
elseif isstring(genes)
    genes = cellstr(genes);
elseif iscategorical(genes)
    genes = cellstr(genes);
end

genes = strtrim(genes(:));

% Convert sample names to a column cell array
if ischar(colnames)
    colnames = cellstr(colnames);
elseif isstring(colnames)
    colnames = cellstr(colnames);
elseif iscategorical(colnames)
    colnames = cellstr(colnames);
end

colnames = colnames(:);

% Check dimensions
if size(data,1) ~= length(genes)
    error('Number of genes does not match rows in data.');
end

if size(data,2) ~= length(colnames)
    error('Number of sample names does not match columns in data.');
end

% Find the selected gene
idxGene = strcmpi(genes, geneName);

if ~any(idxGene)
    error('Gene %s was not found.', geneName);
end

geneIndex = find(idxGene, 1, 'first');

% Expression of selected gene
geneValues = data(geneIndex, :);

% Quartile thresholds
thresholdHigh = prctile(geneValues, 75);
thresholdLow = prctile(geneValues, 25);

% Select samples
selectON = geneValues >= thresholdHigh;
selectOFF = geneValues <= thresholdLow;

% Selected sample names
colnamesON = colnames(selectON);
colnamesOFF = colnames(selectOFF);

% Create gene column
geneTable = table(genes, 'VariableNames', {'Genes'});

% Create HIGH data table
dataTableON = array2table(data(:, selectON));
dataTableON.Properties.VariableNames = colnamesON';

Table_ON = [geneTable, dataTableON];

% Create LOW data table
dataTableOFF = array2table(data(:, selectOFF));
dataTableOFF.Properties.VariableNames = colnamesOFF';

Table_OFF = [geneTable, dataTableOFF];

fprintf('%s HIGH samples: %d\n', geneName, sum(selectON));
fprintf('%s LOW samples: %d\n', geneName, sum(selectOFF));

end
