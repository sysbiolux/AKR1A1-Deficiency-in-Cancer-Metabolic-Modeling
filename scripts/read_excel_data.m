function data = read_excel_data(filename)
%READ_EXCEL_DATA Read an Excel file as a tdfread-style structure.
%
% MATLAB R2019b compatible.

if exist(filename,'file') ~= 2
    error('File not found: %s',filename);
end

[pathstr,name,~] = fileparts(filename);

if isempty(pathstr)
    pathstr = pwd;
end

temp_txt_file = fullfile(pathstr,[name '_temporary.txt']);

opts = detectImportOptions(filename);
T = readtable(filename,opts);

writetable(T, ...
    temp_txt_file, ...
    'FileType','text', ...
    'Delimiter','\t', ...
    'WriteVariableNames',true);

data = tdfread(temp_txt_file,'\t');

if exist(temp_txt_file,'file') == 2
    delete(temp_txt_file);
end

end