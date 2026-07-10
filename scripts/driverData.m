%%
clear; close all; clc;
solverOK=1;
%% Load normalized gene-expression data
% Import normalized count tables from the Zenodo repository data folder.

temp = tdfread(fullfile('data','raw','HuH7', ...
    'Huh7_project_1_normalizedCountsWithAnnotations.txt'));

temp2 = tdfread(fullfile('data','raw','769-p', ...
    '769P_project_1_normalizedCountsWithAnnotations.txt'));


%% Extract sample column names

temp3 = fieldnames(temp);
colnamesH = temp3(8:18);

temp3 = fieldnames(temp2);
colnames7 = temp3(8:18);


%% Identify entries without an annotated gene symbol

isPseudo = find(ismember(cellstr(temp.Symbol),'NA'));
numel(isPseudo)


%% Reorder samples
% Assemble expression matrices using the original sample order.

data_h = [temp.Huh7_siSCR_1, temp.Huh7_siSCR_2, temp.Huh7_siSCR_3, temp.Huh7_siSCR_4];
data_h = [data_h, temp.Huh7_siAKR1A1_1_2, temp.Huh7_siAKR1A1_1_3, temp.Huh7_siAKR1A1_1_4];
data_h = [data_h, temp.Huh7_siAKR1A1_2_1, temp.Huh7_siAKR1A1_2_2, temp.Huh7_siAKR1A1_2_3, temp.Huh7_siAKR1A1_2_4];


data_r = [temp2.x7690x2DP_siSCR_1, temp2.x7690x2DP_siSCR_2, temp2.x7690x2DP_siSCR_3, temp2.x7690x2DP_siSCR_4];
data_r = [data_r, temp2.x7690x2DP_siAKR1A1_1_1, temp2.x7690x2DP_siAKR1A1_1_2, temp2.x7690x2DP_siAKR1A1_1_3];
data_r = [data_r, temp2.x7690x2DP_siAKR1A1_2_1, temp2.x7690x2DP_siAKR1A1_2_2, temp2.x7690x2DP_siAKR1A1_2_3, temp2.x7690x2DP_siAKR1A1_2_4];
%%
%% Output folder
% Create the folder used to save figures if it does not already exist.

if ~exist('images','dir')
    mkdir('images')
end


%% 769-P data

data = data_r;
nansum(data)

% Boxplots of normalized expression values for all samples.
figure
boxplot(data)
saveas(gcf,fullfile('images','r_boxplot_all.png'))
saveas(gcf,fullfile('images','r_boxplot_all.fig'))
close

% Distribution of normalized expression values in the first sample.
figure
histogram(data(:,1),100)
saveas(gcf,fullfile('images','r_gene_distribution_all.png'))
saveas(gcf,fullfile('images','r_gene_distribution_all.fig'))
close

% Kernel density distributions for all samples.
figure
for counter = 1:size(data,2)
    [probability_estimate,xi] = ksdensity(data(:,counter));
    plot(xi,probability_estimate,'-k','LineWidth',1)
    ylabel('Density')
    xlabel('Expression')
    hold on
end
saveas(gcf,fullfile('images','r_cdf_all.png'))
saveas(gcf,fullfile('images','r_cdf_all.fig'))
close


%% HuH7 data

data = data_h;
nansum(data)

% Boxplots of normalized expression values for all samples.
figure
boxplot(data)
saveas(gcf,fullfile('images','h_boxplot_all.png'))
saveas(gcf,fullfile('images','h_boxplot_all.fig'))
close

% Distribution of normalized expression values in the first sample.
figure
histogram(data(:,1),100)
saveas(gcf,fullfile('images','h_gene_distribution_all.png'))
saveas(gcf,fullfile('images','h_gene_distribution_all.fig'))
close

% Kernel density distributions for all samples.
figure
for counter = 1:size(data,2)
    [probability_estimate,xi] = ksdensity(data(:,counter));
    plot(xi,probability_estimate,'-k','LineWidth',1)
    ylabel('Density')
    xlabel('Expression')
    hold on
end
saveas(gcf,fullfile('images','h_cdf_all.png'))
saveas(gcf,fullfile('images','h_cdf_all.fig'))
close
%% PCA
%% PCA: 769-P data
data = data_r;

% Perform PCA using samples as observations.
[~,score,~,~,explained,~] = pca(data');

hold on
plot(score(1:4,1),score(1:4,2),'bo','MarkerFaceColor','b')
plot(score(5:7,1),score(5:7,2),'ro','MarkerFaceColor','r')
plot(score(8:11,1),score(8:11,2),'mo','MarkerFaceColor','m')

title('PCA')
xlabel([num2str(1), ' component: ', num2str(round(explained(1),1))])
ylabel([num2str(2), ' component: ', num2str(round(explained(2),1))])
legend({'SCR','AKR1A1-1','AKR1A1-2'},'location','best')

saveas(gcf,'images/r_pca_score_1.png')
saveas(gcf,'images/r_pca_score_1.fig')
close


%% PCA: HuH7 data
data = data_h;

% Perform PCA using samples as observations.
[~,score,~,~,explained,~] = pca(data');

hold on
plot(score(1:4,1),score(1:4,2),'bo','MarkerFaceColor','b')
plot(score(5:7,1),score(5:7,2),'ro','MarkerFaceColor','r')
plot(score(8:11,1),score(8:11,2),'mo','MarkerFaceColor','m')

title('PCA')
xlabel([num2str(1), ' component: ', num2str(round(explained(1),1))])
ylabel([num2str(2), ' component: ', num2str(round(explained(2),1))])
legend({'SCR','AKR1A1-1','AKR1A1-2'},'location','best')

saveas(gcf,'images/h_pca_score_1.png')
saveas(gcf,'images/h_pca_score_1.fig')
close

%% Remove zero values and visually inspect the rlog distributions

%% 769-P data
sum(data_r==0)

data_non_zeros = data_r;
data_non_zeros(data_non_zeros==0) = NaN;

boxplot(data_non_zeros)
saveas(gcf,'images/r_boxplot.png')
saveas(gcf,'images/r_boxplot.fig')
close

figure
color = 'krcgyb';

for counter = 1:size(data_non_zeros,2)
    hold on
    [probability_estimate,xi] = ksdensity(data_non_zeros(:,counter));

    if counter<=6
        plot(xi,probability_estimate,'-','LineWidth',1,'Color',color(counter));
    else
        plot(xi,probability_estimate,'--','LineWidth',1,'Color',color(counter-6));
    end

    ylabel('Density')
    xlabel('Expression')
end

saveas(gcf,'images/r_cdf.png')
saveas(gcf,'images/r_cdf.fig')
close

histogram(data_non_zeros(:,1),100)
saveas(gcf,'images/r_gene_distribution.png')
saveas(gcf,'images/r_gene_distribution.fig')
close

% Visual inspection of fitted curves for the rlog counts.
for counter = 1:size(data_non_zeros,2)
    fittingcurve_TS(data_non_zeros(:,counter))
end

saveas(gcf,'images/r_fittingcurve_TS.png')
saveas(gcf,'images/r_fittingcurve_TS.fig')
close


%% HuH7 data
sum(data_h==0)

data_non_zeros = data_h;
data_non_zeros(data_non_zeros==0) = NaN;

boxplot(data_non_zeros)
saveas(gcf,'images/h_boxplot.png')
saveas(gcf,'images/h_boxplot.fig')
close

hold on
color = 'krcgyb';

for counter = 1:size(data_non_zeros,2)
    [probability_estimate,xi] = ksdensity(data_non_zeros(:,counter));

    if counter<=6
        plot(xi,probability_estimate,'-','LineWidth',1,'Color',color(counter));
    else
        plot(xi,probability_estimate,'--','LineWidth',1,'Color',color(counter-6));
    end

    ylabel('Density')
    xlabel('Expression')
end

saveas(gcf,'images/h_cdf.png')
saveas(gcf,'images/h_cdf.fig')
close

histogram(data_non_zeros(:,1),100)
saveas(gcf,'images/h_gene_distribution.png')
saveas(gcf,'images/h_gene_distribution.fig')
close

% Visual inspection of fitted curves for the rlog counts.
for counter = 1:size(data_non_zeros,2)
    fittingcurve_TS(data_non_zeros(:,counter))
end

saveas(gcf,'images/h_fittingcurve_TS.png')
saveas(gcf,'images/h_fittingcurve_TS.fig')
close
%%
%% Load gene lengths

% Extract gene IDs and symbols from the imported annotation tables.
genes_id_H = cellstr(temp.Alias);
numel(unique(genes_id_H))

genes_id_7 = cellstr(temp2.Alias);
numel(unique(genes_id_7))

geneSyms_orig_H = cellstr(temp.Symbol);
geneSyms_orig_7 = cellstr(temp2.Symbol);


% Load gene-length information.
[NUM,TXT,RAW] = xlsread(fullfile('data','raw','geneLengths_170423.xlsx'));

lengths_orig = NUM(:,1);
lengths = lengths_orig;

% Check gene-length distribution and missing values.
numel(lengths)
sum(isnan(lengths))
sum(lengths==0)
sum(lengths<100)


% Visualize gene-length distribution.
histogram(log10(lengths),100)
title('gene lengths')
xlabel('log10(gene lengths)')

saveas(gcf,'images/gene_lengths.png')
saveas(gcf,'images/gene_lengths.fig')
close
%%
%% Load gene lengths

% Extract gene IDs and symbols from the imported annotation tables.
genes_id_H = cellstr(temp.Alias);
numel(unique(genes_id_H))

genes_id_7 = cellstr(temp2.Alias);
numel(unique(genes_id_7))

geneSyms_orig_H = cellstr(temp.Symbol);
geneSyms_orig_7 = cellstr(temp2.Symbol);


% Load gene-length information.
[NUM,TXT,RAW] = xlsread(fullfile('data','raw','geneLengths_170423.xlsx'));

lengths_orig = NUM(:,1);
lengths = lengths_orig;

% Check gene-length distribution and missing values.
numel(lengths)
sum(isnan(lengths))
sum(lengths==0)
sum(lengths<100)


% Visualize gene-length distribution.
histogram(log10(lengths),100)
title('gene lengths')
xlabel('log10(gene lengths)')

saveas(gcf,'images/gene_lengths.png')
saveas(gcf,'images/gene_lengths.fig')
close

%% TPM/FPKM scaling

for counter2 = 1:2

    if counter2 == 1
        data = data_h;
        geneSyms = geneSyms_orig_H;
        genes_id = genes_id_H;
    else
        data = data_r;
        geneSyms = geneSyms_orig_7;
        genes_id = genes_id_7;
    end

    % Convert transformed values back to the original scale.
    data = 2.^data - 1;

    % Set pseudo-gene values to zero.
    data(isPseudo,:) = 0;

    lengths = lengths_orig;

    % Remove genes without available gene-length information.
    remove = isnan(lengths);

    if sum(remove) > 0
        data(remove,:) = [];
        lengths(remove) = [];
        genes_id(remove) = [];
        geneSyms(remove) = [];
    end

    % TPM scaling.
    temp = data ./ repmat(lengths,1,size(data,2));
    TPM = 1e6 * temp ./ repmat(nansum(temp),size(data,1),1);

    % FPKM scaling.
    temp = 1e6 * data ./ repmat(nansum(data),size(data,1),1);
    FPKM = temp ./ repmat(lengths,1,size(data,2));

    % Select TPM for the following analyses.
    data3 = TPM;

    nansum(data3)
    sum(data3 == 0)

    data3(data3 == 0) = NaN;
    data3 = log2(data3 + 0);

    % Boxplot of TPM values.
    boxplot(data3)

    if counter2 == 1
        saveas(gcf,'images/boxplot_TPM_H.png')
        saveas(gcf,'images/boxplot_TPM_H.fig')
    else
        saveas(gcf,'images/boxplot_TPM_7.png')
        saveas(gcf,'images/boxplot_TPM_7.fig')
    end

    close

    % Histogram of TPM values for the first sample.
    histogram(data3(:,1),100)

    if counter2 == 1
        saveas(gcf,'images/boxplot_TPM_H.png')
        saveas(gcf,'images/boxplot_TPM_H.fig')
    else
        saveas(gcf,'images/boxplot_TPM_7.png')
        saveas(gcf,'images/boxplot_TPM_7.fig')
    end

    if counter2 == 1
        saveas(gcf,'images/hist_TPM_H.png')
        saveas(gcf,'images/hist_TPM_H.fig')
    else
        saveas(gcf,'images/hist_TPM_7.png')
        saveas(gcf,'images/hist_TPM_7.fig')
    end

    close

    % Density distributions of TPM values.
    figure

    for counter = 1:size(data3,2)
        [probability_estimate,xi] = ksdensity(data3(:,counter));
        plot(xi,probability_estimate,'-k','LineWidth',1);
        ylabel('Density')
        xlabel('expression (TPM/FPKM)')
    end

    if counter2 == 1
        saveas(gcf,'images/density_TPM_H.png')
        saveas(gcf,'images/density_TPM_H.fig')
    else
        saveas(gcf,'images/density_TPM_7.png')
        saveas(gcf,'images/density_TPM_7.fig')
    end

    close

    %% Visual inspection: fit curves to TPM/FPKM values

    for counter = 1:size(data3,2)
        fittingcurve_TS(data3(:,counter))
    end

    saveas(gcf,'images/r_fittingcurve_TS.png')
    saveas(gcf,'images/r_fittingcurve_TS.fig')
    close

    %% Discretize TPM data

    % For data_h:
    % discretized = discretize_FPKM(TPM,colnamesH,1);

    % For data_r:
    discretized = discretize_FPKM(TPM,colnames7,1);

    % Variables to be saved:
    TPM;
    colnamesH;
    colnames7;
    genes_id;

    if counter2 == 1
        save discretizedH TPM geneSyms colnamesH discretized
    else
        save discretized7 TPM geneSyms colnames7 discretized
    end

    %% Statistics for discretized data

    on = sum(discretized == 1);
    disp(on)

    nd = sum(discretized == 0);
    disp(nd)

    off = sum(discretized == -1);
    disp(off)

    tmp = on + nd + off;
    disp(tmp)

end
%%
save driverData
