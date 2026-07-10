function sampleTable = import_TCGA_cancer_type_samples(filename, dataLines)

if nargin < 2
    dataLines = [1, Inf];
end

opts = delimitedTextImportOptions('NumVariables', 2);
opts.DataLines = dataLines;
opts.Delimiter = '\t';

opts.VariableNames = {'SampleID','CancerType'};
opts.VariableTypes = {'char','categorical'};

opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';

opts = setvaropts(opts, 'SampleID', 'WhitespaceRule', 'preserve');
opts = setvaropts(opts, {'SampleID','CancerType'}, 'EmptyFieldRule', 'auto');

sampleTable = readtable(filename, opts);

end