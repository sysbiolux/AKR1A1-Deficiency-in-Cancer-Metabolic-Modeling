function data = read_excel_data(filename)
%READ_EXCEL_DATA Read an Excel file as a tdfread-style structure.

if exist(filename, 'file') ~= 2
    error('File not found: %%s', filename);
end

[pathstr, name, ~] = fileparts(filename);

temp_txt_file = fullfile(pathstr, [name '_temporary.txt']);

opts = detectImportOptions(filename);
T = readtable(filename, opts);

writetable(T, ...
    temp_txt_file, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'WriteVariableNames', true);

data = tdfread(temp_txt_file, '\t');

if exist(temp_txt_file, 'file') == 2
    delete(temp_txt_file);
end

end
