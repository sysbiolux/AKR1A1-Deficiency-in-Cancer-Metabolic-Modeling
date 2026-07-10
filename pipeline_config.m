function cfg = pipeline_config()

rootDir = fileparts(mfilename('fullpath'));

cfg.rootDir = rootDir;
cfg.dataDir = fullfile(rootDir, 'data');

cfg.rawTCGADir = fullfile(cfg.dataDir, 'raw_TCGA');
cfg.rawCellDir = fullfile(cfg.dataDir, 'raw_cell_lines');
cfg.cancerExpressionDir = fullfile(cfg.dataDir, 'cancer_expression_TCGA');
cfg.highLowDir = fullfile(cfg.dataDir, 'high_low');

cfg.modelDir = fullfile(cfg.dataDir, 'model');
cfg.mediumDir = fullfile(cfg.dataDir, 'medium');
cfg.drugDir = fullfile(cfg.dataDir, 'drug');

cfg.resultsDir = fullfile(rootDir, 'results');
cfg.preparedDir = fullfile(cfg.resultsDir, '01_prepared_data');
cfg.modelsDir = fullfile(cfg.resultsDir, '02_models');
cfg.analysisDir = fullfile(cfg.resultsDir, '03_analysis');
cfg.figuresDir = fullfile(cfg.resultsDir, 'figures');

% Choose: 'TCGA', 'CELL_LINES', or 'ALL'
cfg.runMode = 'TCGA';

% TCGA GEO GSE62944 files
cfg.tcga.sampleAnnotationFile = 'GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt';
cfg.tcga.tumorExpressionFile = 'GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt';
cfg.tcga.normalExpressionFile = 'GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FPKM.txt';

cfg.cancerTypes = {'KICH','KIRC','KIRP','LIHC'};

cfg.highLow.gene = 'AKR1A1';
cfg.highLow.highPercentile = 75;
cfg.highLow.lowPercentile = 25;

% Cell-line GEO folders
cfg.cellLineDir.GSE310784 = fullfile(cfg.rawCellDir, 'GSE310784', '769-p');
cfg.cellLineDir.GSE310828 = fullfile(cfg.rawCellDir, 'GSE310828', 'HuH7');

cfg.cellLines.H.name = 'Huh7';
cfg.cellLines.H.expression = fullfile(cfg.cellLineDir.GSE310828, 'Huh7_project_1_normalizedCountsWithAnnotations.txt');
cfg.cellLines.H.medium = 'DMEM';
cfg.cellLines.H.groupNames = {'ctrl','sc1','sc2','sc12'};
cfg.cellLines.H.groupColumns = {[1 2 3 11], 4:6, 7:10, 4:10};

cfg.cellLines.Seven.name = '769p';
cfg.cellLines.Seven.expression = fullfile(cfg.cellLineDir.GSE310784, '769P_project_1_normalizedCountsWithAnnotations.txt');
cfg.cellLines.Seven.medium = 'RPMI';
cfg.cellLines.Seven.groupNames = {'ctrl','sc1','sc2','sc12'};
cfg.cellLines.Seven.groupColumns = {[1 2 3 4], [5 6 7], [8 9 10 11], 5:11};

% Model support files
cfg.modelFile = fullfile(cfg.modelDir, 'consistent_model.mat');
cfg.dicoFile = fullfile(cfg.modelDir, 'dico2columns.mat');

cfg.mediumFiles.DMEM = fullfile(cfg.mediumDir, 'medium_DMEM_H.xlsx');
cfg.mediumFiles.RPMI = fullfile(cfg.mediumDir, 'medium_RPMI_7.xlsx');

cfg.drugFile = fullfile(cfg.drugDir, 'GeneDrugRelationsUpdate.mat');

% FASTCORMICS settings
cfg.fastcormics.already_mapped_tag = 0;
cfg.fastcormics.consensus_proportion = 0.9;
cfg.fastcormics.epsilon = 1e-4;
cfg.fastcormics.biomass_rxn = 'biomass_reaction';
cfg.fastcormics.func = {'DM_atp_c_', 'biomass_reaction'};

cfg.fastcormics.unpenalizedSystems = { ...
    'Transport, endoplasmic reticular'; ...
    'Transport, extracellular'; ...
    'Transport, golgi apparatus'; ...
    'Transport, mitochondrial'; ...
    'Transport, peroxisomal'; ...
    'Transport, lysosomal'; ...
    'Transport, nuclear'};

cfg.mediumByCancer.KICH = 'RPMI';
cfg.mediumByCancer.KIRC = 'RPMI';
cfg.mediumByCancer.KIRP = 'RPMI';
cfg.mediumByCancer.LIHC = 'DMEM';

folders = {cfg.dataDir, cfg.rawTCGADir, cfg.rawCellDir, ...
    cfg.cancerExpressionDir, cfg.highLowDir, cfg.modelDir, ...
    cfg.mediumDir, cfg.drugDir, cfg.resultsDir, cfg.preparedDir, ...
    cfg.modelsDir, cfg.analysisDir, cfg.figuresDir};

for i = 1:numel(folders)
    if ~exist(folders{i}, 'dir')
        mkdir(folders{i});
    end
end

end