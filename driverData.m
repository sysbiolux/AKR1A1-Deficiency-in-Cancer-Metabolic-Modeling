%% driverData_highlow_clean.m
% Data driver for both:
% 1) Cell lines: Huh7 and 769-p
% 2) TCGA AKR1A1 ON/OFF: KICH, KIRC, KIRP, LIHC
%
% MATLAB R2019b compatible.
% Keeps original logic and uses discretize_FPKM.

clear; close all; clc;
solverOK = 1;

cfg = pipeline_config;

if ~exist(fullfile(cfg.rootDir,'images'),'dir')
    mkdir(fullfile(cfg.rootDir,'images'));
end

%% ============================================================
% TCGA HIGH / LOW DATA
%% ============================================================

if strcmpi(cfg.runMode,'TCGA') || strcmpi(cfg.runMode,'ALL')

    files_off = {'Table_OFF_KICH','Table_OFF_KIRC','Table_OFF_KIRP','Table_OFF_LIHC'};
    files_on  = {'Table_ON_KICH', 'Table_ON_KIRC', 'Table_ON_KIRP', 'Table_ON_LIHC'};
    cancers   = {'KICH','KIRC','KIRP','LIHC'};

    for i = 1:length(cancers)

        cancer = cancers{i};
        disp(['Preparing TCGA ' cancer]);

        file_off = fullfile(cfg.highLowDir,[files_off{i} '.xlsx']);
        file_on  = fullfile(cfg.highLowDir,[files_on{i} '.xlsx']);

        temp_off = read_excel_data(file_off);
        temp_on  = read_excel_data(file_on);

        temp = temp_off;
        temp2 = temp_on;

        temp3 = fieldnames(temp);
        colnamesOFF = temp3(2:end);

        temp3 = fieldnames(temp2);
        colnamesON = temp3(2:end);

        data_valuesoff = [];
        for k = 1:length(colnamesOFF)
            data_valuesoff = [data_valuesoff, double(temp.(colnamesOFF{k}))];
        end

        data_valueson = [];
        for k = 1:length(colnamesON)
            data_valueson = [data_valueson, double(temp2.(colnamesON{k}))];
        end

        symbolFieldOff = get_symbol_field(temp);
        symbolFieldOn  = get_symbol_field(temp2);

        geneSymsoff = cellstr(temp.(symbolFieldOff));
        geneSymson  = cellstr(temp2.(symbolFieldOn));

        isPseudoOff = find(ismember(geneSymsoff,'NA'));
        isPseudoOn  = find(ismember(geneSymson,'NA'));

        disp(['Pseudo genes OFF: ' num2str(numel(isPseudoOff))]);
        disp(['Pseudo genes ON: ' num2str(numel(isPseudoOn))]);

        discretizedon  = discretize_FPKM(data_valueson,  colnamesON,  1);
        discretizedoff = discretize_FPKM(data_valuesoff, colnamesOFF, 1);

        save(fullfile(cfg.preparedDir,['discretizedON_' cancer '.mat']), ...
            'discretizedon','geneSymson','colnamesON');

        save(fullfile(cfg.preparedDir,['discretizedOFF_' cancer '.mat']), ...
            'discretizedoff','geneSymsoff','colnamesOFF');

        on = sum(discretizedon == 1);
        nd = sum(discretizedon == 0);
        off = sum(discretizedon == -1);
        disp(['ON stats for ' cancer]);
        disp(on); disp(nd); disp(off); disp(on+nd+off);

        on = sum(discretizedoff == 1);
        nd = sum(discretizedoff == 0);
        off = sum(discretizedoff == -1);
        disp(['OFF stats for ' cancer]);
        disp(on); disp(nd); disp(off); disp(on+nd+off);

        save(fullfile(cfg.preparedDir,['driverDataOFFON_' cancer '.mat']), ...
            'data_valueson','data_valuesoff', ...
            'geneSymson','geneSymsoff', ...
            'colnamesON','colnamesOFF', ...
            'discretizedon','discretizedoff');
    end
end

%% ============================================================
% CELL LINE DATA: Huh7 and 769-p
%% ============================================================

if strcmpi(cfg.runMode,'CELL_LINES') || strcmpi(cfg.runMode,'ALL')

    huh7File = fullfile(cfg.cellLineDir.Huh7, ...
        'Huh7_project_1_normalizedCountsWithAnnotations.txt');

    p769File = fullfile(cfg.cellLineDir.p769, ...
        '769P_project_1_normalizedCountsWithAnnotations.txt');

    temp  = tdfread(huh7File);
    temp2 = tdfread(p769File);

    temp3 = fieldnames(temp);
    colnamesH = temp3(8:18);

    temp3 = fieldnames(temp2);
    colnames7 = temp3(8:18);

    isPseudoH = find(ismember(cellstr(temp.Symbol),'NA'));
    isPseudo7 = find(ismember(cellstr(temp2.Symbol),'NA'));

    disp(['Pseudo genes Huh7: ' num2str(numel(isPseudoH))]);
    disp(['Pseudo genes 769-p: ' num2str(numel(isPseudo7))]);

    %% Original sample reordering

    data_h = [temp.Huh7_siSCR_1, temp.Huh7_siSCR_2, temp.Huh7_siSCR_3, temp.Huh7_siSCR_4];
    data_h = [data_h, temp.Huh7_siAKR1A1_1_2, temp.Huh7_siAKR1A1_1_3, temp.Huh7_siAKR1A1_1_4];
    data_h = [data_h, temp.Huh7_siAKR1A1_2_1, temp.Huh7_siAKR1A1_2_2, temp.Huh7_siAKR1A1_2_3, temp.Huh7_siAKR1A1_2_4];

    data_r = [temp2.x7690x2DP_siSCR_1, temp2.x7690x2DP_siSCR_2, temp2.x7690x2DP_siSCR_3, temp2.x7690x2DP_siSCR_4];
    data_r = [data_r, temp2.x7690x2DP_siAKR1A1_1_1, temp2.x7690x2DP_siAKR1A1_1_2, temp2.x7690x2DP_siAKR1A1_1_3];
    data_r = [data_r, temp2.x7690x2DP_siAKR1A1_2_1, temp2.x7690x2DP_siAKR1A1_2_2, temp2.x7690x2DP_siAKR1A1_2_3, temp2.x7690x2DP_siAKR1A1_2_4];

    genes_id_H = cellstr(temp.Alias);
    genes_id_7 = cellstr(temp2.Alias);

    geneSyms_orig_H = cellstr(temp.Symbol);
    geneSyms_orig_7 = cellstr(temp2.Symbol);

    geneLengthFile = fullfile(cfg.dataDir,'geneLengths_170423.xlsx');

    if exist(geneLengthFile,'file') ~= 2
        geneLengthFile = fullfile(cfg.rootDir,'geneLengths_170423.xlsx');
    end

    [NUM,TXT,RAW] = xlsread(geneLengthFile);
    lengths_orig = NUM(:,1);

    for counter2 = 1:2

        if counter2 == 1
            data = data_h;
            geneSyms = geneSyms_orig_H;
            genes_id = genes_id_H;
            colnamesCurrent = colnamesH;
            isPseudo = isPseudoH;
            saveName = 'discretizedH.mat';
            label = 'H';
        else
            data = data_r;
            geneSyms = geneSyms_orig_7;
            genes_id = genes_id_7;
            colnamesCurrent = colnames7;
            isPseudo = isPseudo7;
            saveName = 'discretized7.mat';
            label = '7';
        end

        data = 2.^data - 1;
        data(isPseudo,:) = 0;

        lengths = lengths_orig;

        remove = isnan(lengths);

        if sum(remove) > 0
            data(remove,:) = [];
            lengths(remove) = [];
            genes_id(remove) = [];
            geneSyms(remove) = [];
        end

        tempTPM = data ./ repmat(lengths,1,size(data,2));
        TPM = 1e+6 * tempTPM ./ repmat(nansum(tempTPM),size(data,1),1);

        tempFPKM = 1e+6 * data ./ repmat(nansum(data),size(data,1),1);
        FPKM = tempFPKM ./ repmat(lengths,1,size(data,2));

        data3 = TPM;
        data3(data3 == 0) = NaN;
        data3 = log2(data3 + 0);

        discretized = discretize_FPKM(TPM, colnamesCurrent, 1);

        if strcmp(label,'H')
            colnamesH = colnamesCurrent;
            save(fullfile(cfg.preparedDir,saveName), ...
                'TPM','FPKM','geneSyms','genes_id','colnamesH','discretized');
        else
            colnames7 = colnamesCurrent;
            save(fullfile(cfg.preparedDir,saveName), ...
                'TPM','FPKM','geneSyms','genes_id','colnames7','discretized');
        end

        on = sum(discretized == 1);
        nd = sum(discretized == 0);
        off = sum(discretized == -1);

        disp(['Discretization stats ' label]);
        disp(on); disp(nd); disp(off); disp(on+nd+off);
    end

    save(fullfile(cfg.preparedDir,'driverData_cell_lines.mat'), ...
        'data_h','data_r','colnamesH','colnames7', ...
        'geneSyms_orig_H','geneSyms_orig_7');
end

save(fullfile(cfg.preparedDir,'driverData_highlow_clean_workspace.mat'));

disp('driverData_highlow_clean finished');


function symbolField = get_symbol_field(S)

fields = fieldnames(S);

if ismember('Symbol',fields)
    symbolField = 'Symbol';
elseif ismember('Var1',fields)
    symbolField = 'Var1';
elseif ismember('Genes',fields)
    symbolField = 'Genes';
else
    symbolField = fields{1};
end

end