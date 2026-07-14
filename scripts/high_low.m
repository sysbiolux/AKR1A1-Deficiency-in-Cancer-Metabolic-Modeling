function [Table_ON, Table_OFF, colnamesON, colnamesOFF] = ...
    high_low(data, genes, colnames, geneName)
%HIGH_LOW Divide samples into HIGH and LOW groups for a selected gene.
%
% HIGH samples:
%   expression >= 75th percentile
%
% LOW samples:
%   expression <= 25th percentile
%
% Inputs
% ------
% data:
%   Numeric genes-by-samples expression matrix.
%
% genes:
%   Gene symbols corresponding to rows of data.
%
% colnames:
%   Sample names corresponding to columns of data.
%
% geneName:
%   Gene used for splitting, for example 'AKR1A1'.
%
% Outputs
% -------
% Table_ON:
%   Gene column followed by HIGH-expression samples.
%
% Table_OFF:
%   Gene column followed by LOW-expression samples.
%
% colnamesON:
%   Names of HIGH-expression samples.
%
% colnamesOFF:
%   Names of LOW-expression samples.
%
% MATLAB R2019b compatible.

    %% Validate expression matrix

    if ~isnumeric(data)
        error('The expression data must be a numeric matrix.');
    end

    if isempty(data)
        error('The expression matrix is empty.');
    end

    %% Convert genes into a column cell array

    if ischar(genes)
        genes = cellstr(genes);
    elseif isstring(genes)
        genes = cellstr(genes);
    elseif iscategorical(genes)
        genes = cellstr(genes);
    end

    genes = strtrim(genes(:));

    %% Convert sample names into a column cell array

    if ischar(colnames)
        colnames = cellstr(colnames);
    elseif isstring(colnames)
        colnames = cellstr(colnames);
    elseif iscategorical(colnames)
        colnames = cellstr(colnames);
    end

    colnames = strtrim(colnames(:));

    %% Check dimensions

    if size(data,1) ~= numel(genes)
        error( ...
            ['Number of genes (%d) does not match the number ' ...
             'of rows in data (%d).'], ...
            numel(genes), ...
            size(data,1));
    end

    if size(data,2) ~= numel(colnames)
        error( ...
            ['Number of sample names (%d) does not match the ' ...
             'number of columns in data (%d).'], ...
            numel(colnames), ...
            size(data,2));
    end

    %% Find the selected gene

    idxGene = strcmpi(genes,geneName);
    numberOfMatches = sum(idxGene);

    if numberOfMatches == 0
        error('Gene %s was not found.',geneName);
    end

    if numberOfMatches > 1
        warning( ...
            ['Gene %s occurs %d times. The first occurrence ' ...
             'will be used for the HIGH/LOW split.'], ...
            geneName, ...
            numberOfMatches);
    end

    geneIndex = find(idxGene,1,'first');

    %% Extract expression of the selected gene

    geneValues = double(data(geneIndex,:));

    validValues = ~isnan(geneValues) & isfinite(geneValues);

    if sum(validValues) < 4
        error( ...
            ['Gene %s has fewer than four finite expression ' ...
             'values. Quartile splitting cannot be performed.'], ...
            geneName);
    end

    %% Calculate quartile thresholds

    thresholdHigh = prctile(geneValues(validValues),75);
    thresholdLow  = prctile(geneValues(validValues),25);

    fprintf('%s 25th percentile: %g\n',geneName,thresholdLow);
    fprintf('%s 75th percentile: %g\n',geneName,thresholdHigh);

    %% Select samples

    selectON = validValues & geneValues >= thresholdHigh;
    selectOFF = validValues & geneValues <= thresholdLow;

    %% Check whether groups overlap

    overlappingSamples = selectON & selectOFF;

    if any(overlappingSamples)
        error( ...
            ['The HIGH and LOW groups overlap for gene %s. ' ...
             'This usually means the 25th and 75th percentile ' ...
             'thresholds are identical.'], ...
            geneName);
    end

    if ~any(selectON)
        error('No HIGH samples were selected for gene %s.',geneName);
    end

    if ~any(selectOFF)
        error('No LOW samples were selected for gene %s.',geneName);
    end

    %% Extract original selected sample names

    colnamesON_original = colnames(selectON);
    colnamesOFF_original = colnames(selectOFF);

    %% Make names valid for MATLAB tables

    colnamesON = matlab.lang.makeValidName(colnamesON_original);
    colnamesOFF = matlab.lang.makeValidName(colnamesOFF_original);

    colnamesON = matlab.lang.makeUniqueStrings(colnamesON);
    colnamesOFF = matlab.lang.makeUniqueStrings(colnamesOFF);

    colnamesON = colnamesON(:);
    colnamesOFF = colnamesOFF(:);

    if ~isequal(colnamesON,colnamesON_original)
        warning('Some HIGH sample names were adjusted for MATLAB tables.');
    end

    if ~isequal(colnamesOFF,colnamesOFF_original)
        warning('Some LOW sample names were adjusted for MATLAB tables.');
    end

    %% Create the gene column

    geneTable = table(genes,'VariableNames',{'Genes'});

    %% Create HIGH-expression table

    dataTableON = array2table(data(:,selectON));
    dataTableON.Properties.VariableNames = colnamesON';

    Table_ON = [geneTable,dataTableON];

    %% Create LOW-expression table

    dataTableOFF = array2table(data(:,selectOFF));
    dataTableOFF.Properties.VariableNames = colnamesOFF';

    Table_OFF = [geneTable,dataTableOFF];

    %% Report selected groups

    fprintf('%s HIGH samples: %d\n',geneName,sum(selectON));
    fprintf('%s LOW samples: %d\n',geneName,sum(selectOFF));

end