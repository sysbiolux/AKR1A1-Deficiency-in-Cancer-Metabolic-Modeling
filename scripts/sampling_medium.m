%% sampling_medium.m
% Unified model-sampling script for:
%
%   1. Standard medium models
%   2. ON/OFF medium models
%
%
% MATLAB compatibility: MATLAB R2019.

%% Select the dataset to sample

% Use exactly one of these values:
datasetType = 'standard';
% datasetType = 'onoff';


%% Relative project paths
% The launcher must start MATLAB from the scripts folder.
%
% Expected structure:
%
%   scripts/
%       sampling_medium.m
%       parsave_custom.m
%
%   data/
%       standard/
%           models_medium.mat
%           sampling/
%
%       onoff/
%           models_mediumonoff.mat
%           sampling/

dataRoot = '../data';


%% Select the input and output folders

if strcmp(datasetType,'standard')

    inputFolder = fullfile(dataRoot,'standard');

    inputFile = fullfile( ...
        inputFolder, ...
        'models_medium.mat');

    outputFolder = fullfile( ...
        inputFolder, ...
        'sampling');


elseif strcmp(datasetType,'onoff')

    inputFolder = fullfile(dataRoot,'onoff');

    inputFile = fullfile( ...
        inputFolder, ...
        'models_mediumonoff.mat');

    outputFolder = fullfile( ...
        inputFolder, ...
        'sampling');


else

    error(['Unknown datasetType: ' datasetType ...
        '. Use ''standard'' or ''onoff''.']);

end


%% Confirm that the input file exists

if ~exist(inputFile,'file')

    error(['Input file not found: ' inputFile]);

end


%% Create the output folder when necessary

if ~exist(outputFolder,'dir')

    mkdir(outputFolder);

end


%% Load the selected model collection

if strcmp(datasetType,'standard')

    load(inputFile,'models_medium');

    % Preserve the original model order.
    data = models_medium;


elseif strcmp(datasetType,'onoff')

    load(inputFile,'models_mediumonoff');

    % Preserve the original model order.
    data = models_mediumonoff;

end


%% Configure cluster characteristics
% Preserve the original number of workers and model count.

maxWorkers = 8;

cluster = parcluster();

cluster.NumWorkers = maxWorkers;


%% Sample the eight models in parallel

parfor i = 1:maxWorkers

    %% Initialize CPLEX and COBRA Toolbox on every worker
    % These paths correspond to the HPC software installation and are kept
    % as in the original sampling scripts.

    addpath(genpath( ...
        '/opt/apps/resif/aion/2020a/epyc/software/CPLEX'));

    addpath(genpath( ...
        '../cobratoolbox'));

    changeCobraSolver('ibm_cplex');


    %% Preserve the original initialization delay

    pause(30);


    %% Retrieve the current model
    % The model ordering is inherited directly from the loaded collection:
    %
    %   data(1).models
    %   ...
    %   data(8).models
    %
    % No model reordering is performed.

    this_model = data(i).models;


    disp(strcat( ...
        'Hi, I am worker #', ...
        num2str(i), ...
        ', and I am happy to do your job! Below is my model.'));

    disp(this_model);


    %% Sampling options
    % These values are unchanged from the original standard and ON/OFF
    % sampling scripts.

    options = [];

    options.nPointsReturned = 1000;

    options.nFiles = 10;

    options.maxTime = 36000;

    options.nWarmupPoints = ...
        2 * size(this_model.S,2);

    options.nStepsPerPoint = ...
        size(this_model.S,2);


    %% Perform ACHR sampling
    % Preserve the same model name passed to sampleCbModel.

    name = strcat( ...
        'model_', ...
        num2str(i));


    [modelSampling,samples] = sampleCbModel( ...
        this_model, ...
        name, ...
        'ACHR', ...
        options);


    %% Save the standard sampling result

    if strcmp(datasetType,'standard')

        % Preserve the original standard structure name and field order.

        samplingResults_medium(i).samples = samples;

        samplingResults_medium(i).modelSampling = ...
            modelSampling;


        outputFile = fullfile( ...
            outputFolder, ...
            strcat( ...
                'samplingResults_medium_', ...
                num2str(options.nPointsReturned), ...
                '_model_', ...
                num2str(i), ...
                '.mat'));


        parsave_custom( ...
            outputFile, ...
            samplingResults_medium(i));


    %% Save the ON/OFF sampling result

    elseif strcmp(datasetType,'onoff')

        % Preserve the original ON/OFF structure name and field order.

        samplingResults_mediumonoff(i).samples = samples;

        samplingResults_mediumonoff(i).modelSampling = ...
            modelSampling;


        outputFile = fullfile( ...
            outputFolder, ...
            strcat( ...
                'samplingResults_mediumonoff_', ...
                num2str(options.nPointsReturned), ...
                '_model_', ...
                num2str(i), ...
                '.mat'));


        parsave_custom( ...
            outputFile, ...
            samplingResults_mediumonoff(i));

    end


    disp(['Saved sampling result: ' outputFile]);

end


%% End the non-interactive MATLAB session

quit
