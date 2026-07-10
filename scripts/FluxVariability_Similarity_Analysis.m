%% ========================================================================
%  Similarity Based on Flux Variability Analysis
%  ========================================================================

%% Flux Variability Analysis - 769-P models

[minFlux_7contr, maxFlux_7contr] = fluxVariability(model_7_ctrl);
[minFlux_7sc1,   maxFlux_7sc1]   = fluxVariability(model_7_sc1);
[minFlux_7sc12,  maxFlux_7sc12]  = fluxVariability(model_7_sc12);
[minFlux_7sc2,   maxFlux_7sc2]   = fluxVariability(model_7_sc2);


%% Store 769-P FVA results in tables

T_fluxes_7ctrl = table( ...
    model_7_ctrl.rxns, ...
    model_7_ctrl.rxnNames, ...
    minFlux_7contr, ...
    maxFlux_7contr);

T_fluxes_7sc1 = table( ...
    model_7_sc1.rxns, ...
    model_7_sc1.rxnNames, ...
    minFlux_7sc1, ...
    maxFlux_7sc1);

T_fluxes_7sc2 = table( ...
    model_7_sc2.rxns, ...
    model_7_sc2.rxnNames, ...
    maxFlux_7sc2);

T_fluxes_7sc12 = table( ...
    model_7_sc12.rxns, ...
    model_7_sc12.rxnNames, ...
    maxFlux_7sc12);


%% Flux Variability Analysis - Huh7 models

[minFlux_Hcontr, maxFlux_Hcontr] = fluxVariability(model_H_ctrl);
[minFlux_Hsc1,   maxFlux_Hsc1]   = fluxVariability(model_H_sc1);
[minFlux_Hsc2,   maxFlux_Hsc2]   = fluxVariability(model_H_sc2);
[minFlux_Hsc12,  maxFlux_Hsc12]  = fluxVariability(model_H_sc12);


%% Store Huh7 FVA results in tables

T_fluxes_Hctrl = table( ...
    model_H_ctrl.rxns, ...
    model_H_ctrl.rxnNames, ...
    minFlux_Hcontr, ...
    maxFlux_Hcontr);

T_fluxes_Hsc1 = table( ...
    model_H_sc1.rxns, ...
    model_H_sc1.rxnNames, ...
    minFlux_Hsc1, ...
    maxFlux_Hsc1);

T_fluxes_Hsc2 = table( ...
    model_H_sc2.rxns, ...
    model_H_sc2.rxnNames, ...
    maxFlux_Hsc2);

T_fluxes_Hsc12 = table( ...
    model_H_sc12.rxns, ...
    model_H_sc12.rxnNames, ...
    maxFlux_Hsc12);


save analysis_part1


%% ========================================================================
%  Map 769-P FVA Results to the Consistent Model
%  ========================================================================

minFlux_7_keep = zeros(numel(consistent_model.rxns), 4);
maxFlux_7_keep = zeros(numel(consistent_model.rxns), 4);


% Control model

minFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_ctrl.rxns), 1) = ...
    minFlux_7contr;

maxFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_ctrl.rxns), 1) = ...
    maxFlux_7contr;


% SC1 model

minFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc1.rxns), 2) = ...
    minFlux_7sc1;

maxFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc1.rxns), 2) = ...
    maxFlux_7sc1;


% SC2 model

minFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc2.rxns), 3) = ...
    minFlux_7sc2;

maxFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc2.rxns), 3) = ...
    maxFlux_7sc2;


% SC12 model

minFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc12.rxns), 4) = ...
    minFlux_7sc12;

maxFlux_7_keep( ...
    ismember(consistent_model.rxns, model_7_sc12.rxns), 4) = ...
    maxFlux_7sc12;


%% Whole-model FVA similarity - 769-P

res_keep = zeros( ...
    size(minFlux_7_keep, 2), ...
    size(minFlux_7_keep, 2));

for i = 1:size(minFlux_7_keep, 2)

    for j = 1:size(minFlux_7_keep, 2)

        v1mins = minFlux_7_keep(:, i);
        v2mins = minFlux_7_keep(:, j);

        v1maxs = maxFlux_7_keep(:, i);
        v2maxs = maxFlux_7_keep(:, j);

        res = FVA_similarity_Thomas( ...
            v1mins, ...
            v1maxs, ...
            v2mins, ...
            v2maxs);

        res_keep(i, j) = res;

    end

end


%% Subsystem-level FVA similarity - 769-P

k = 0;
res2_keep = zeros(1, 3);

% disp('contol-s1', 'contol-s2', 's1-s2')

for i = 1:size(minFlux_7_keep, 2)

    for j = 1:size(minFlux_7_keep, 2)

        if i > j

            k = k + 1;

            for counter = 1:numel(uniSys)

                match = ismember( ...
                    consistent_model.subSystems, ...
                    uniSys(counter));

                v1mins = minFlux_7_keep(match, i);
                v2mins = minFlux_7_keep(match, j);

                v1maxs = maxFlux_7_keep(match, i);
                v2maxs = maxFlux_7_keep(match, j);

                res = FVA_similarity_Thomas( ...
                    v1mins, ...
                    v1maxs, ...
                    v2mins, ...
                    v2maxs);

                res2_keep(counter, k) = res;

            end

        end

    end

end


T_FVA_7 = table(uniSys, res2_keep);

save('T_FVA_7.mat', 'T_FVA_7');


%% ========================================================================
%  Map Huh7 FVA Results to the Consistent Model
%  ========================================================================

minFlux_H_keep = zeros(numel(consistent_model.rxns), 4);
maxFlux_H_keep = zeros(numel(consistent_model.rxns), 4);


% Control model

minFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_ctrl.rxns), 1) = ...
    minFlux_Hcontr;

maxFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_ctrl.rxns), 1) = ...
    maxFlux_Hcontr;


% SC1 model

minFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc1.rxns), 2) = ...
    minFlux_Hsc1;

maxFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc1.rxns), 2) = ...
    maxFlux_Hsc1;


% SC2 model

minFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc2.rxns), 3) = ...
    minFlux_Hsc2;

maxFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc2.rxns), 3) = ...
    maxFlux_Hsc2;


% SC12 model

minFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc12.rxns), 4) = ...
    minFlux_Hsc12;

maxFlux_H_keep( ...
    ismember(consistent_model.rxns, model_H_sc12.rxns), 4) = ...
    maxFlux_Hsc12;


%% Whole-model FVA similarity - Huh7

res_keep = zeros( ...
    size(minFlux_H_keep, 2), ...
    size(minFlux_H_keep, 2));

for i = 1:size(minFlux_H_keep, 2)

    for j = 1:size(minFlux_H_keep, 2)

        v1mins = minFlux_H_keep(:, i);
        v2mins = minFlux_H_keep(:, j);

        v1maxs = maxFlux_H_keep(:, i);
        v2maxs = maxFlux_H_keep(:, j);

        res = FVA_similarity_Thomas( ...
            v1mins, ...
            v1maxs, ...
            v2mins, ...
            v2maxs);

        res_keep(i, j) = res;

    end

end


%% Subsystem-level FVA similarity - Huh7

k = 0;
res2_keep = zeros(1, 3);

% disp('contol-s1'; 'contol-s2'; 's1-s2');

uniSys = unique(consistent_model.subSystems);

for i = 1:size(minFlux_H_keep, 2)

    for j = 1:size(minFlux_H_keep, 2)

        if i > j

            k = k + 1;

            for counter = 1:numel(uniSys)

                match = ismember( ...
                    consistent_model.subSystems, ...
                    uniSys(counter));

                v1mins = minFlux_H_keep(match, i);
                v2mins = minFlux_H_keep(match, j);

                v1maxs = maxFlux_H_keep(match, i);
                v2maxs = maxFlux_H_keep(match, j);

                res = FVA_similarity_Thomas( ...
                    v1mins, ...
                    v1maxs, ...
                    v2mins, ...
                    v2maxs);

                res2_keep(counter, k) = res;

            end

        end

    end

end


T_FVA_H = table(uniSys, res2_keep);

save('T_FVA_H.mat', 'T_FVA_H');


%% ========================================================================
%  FVA Similarity Heatmaps
%  ========================================================================

T_FVA_7_decrese = readtable( ...
    'FVA_7_decrese_sc12.xlsx', ...
    'PreserveVariableNames', true);

T_FVA_H_decrese = readtable( ...
    'FVA_H_decrese_sc12.xlsx', ...
    'PreserveVariableNames', true);


%% Extract similarity matrices

matSignalT7 = T_FVA_7_decrese(1:30, 2:4);
matSignalTH = T_FVA_H_decrese(1:30, 2:4);

matSignal = [matSignalT7, matSignalTH];
matSignal = table2array(matSignal);

custom_colormap = redbluecmap(101);

matSignalT7 = table2array(matSignalT7);
matSignalTH = table2array(matSignalTH);


%% Extract subsystem names

unisubSys_name_7 = T_FVA_7_decrese.unisubSys(1:30);
unisubSys_name_H = T_FVA_H_decrese.unisubSys(1:30);


%% Define output names and model labels

name_7_RCC = strcat('Similarity_FVA_all_subsystems_RCC');
name_H_HCC = strcat('Similarity_FVA_all_subsystems_HCC');


modelNames_7_RCC = { ...
    'RCC\_ctrl\_sc1', ...
    'RCC\_ctrl\_sc2', ...
    'RCC\_ctrl\_sc12'};

modelNames_H_HCC = { ...
    'HCC\_ctrl\_sc1', ...
    'HCC\_ctrl\_sc2', ...
    'HCC\_ctrl\_sc12'};


modelNames_7_decrease = { ...
    'model_769-P_sc1', ...
    'model_769-P_sc2', ...
    'model_769-P_sc12'};

modelNames_H_decrease = { ...
    'model_Huh7_sc1', ...
    'model_Huh7_sc2', ...
    'model_Huh7_sc12'};


modelNames_7_decrease2 = { ...
    'model_7_sc1', ...
    'model_7_sc2', ...
    'model_7_sc12'};

modelNames_H_decrease2 = { ...
    'model_H_sc1', ...
    'model_H_sc2', ...
    'model_H_sc12'};


feature astheightlimit 2000


%% Define heatmap colormap

altcolor = [ ...
    255 255 255;
    255 204 204;
    255 153 153;
    255 102 102;
    255 51 51;
    255 0 0;
    204 0 0;
    152 0 0;
    102 0 0;
    51 0 0] / 255;


%% FVA similarity heatmap - RCC / 769-P

cgo_J_7 = clustergram( ...
    matSignalT7, ...
    'RowLabels', unisubSys_name_7, ...
    'ColumnLabels', modelNames_7_decrease, ...
    'ColumnLabelsRotate', 360, ...
    'Cluster', 'all', ...
    'symmetric', false, ...
    'Colormap', altcolor);


addTitle( ...
    cgo_J_7, ...
    'Similarity based on Flux Variability Analysis RCC');


drawnow;


cgf = plot(cgo_J_7);


colorbar(cgf, 'eastoutside');


saveas(cgf, 'FVA_similarity_heatmap_7_.png');
saveas(cgf, 'FVA_similarity_heatmap_7_.pdf');


%% FVA similarity heatmap - HCC / Huh7

cgo_J_H = clustergram( ...
    matSignalTH(:, :), ...
    'RowLabels', unisubSys_name_H, ...
    'ColumnLabels', modelNames_H_decrease, ...
    'ColumnLabelsRotate', 360, ...
    'Cluster', 'all', ...
    'symmetric', 'False', ...
    'Colormap', altcolor);


addTitle( ...
    cgo_J_H, ...
    'Similarity based on Flux Variability Analysis HCC');


drawnow;


cgf = plot(cgo_J_H);


colorbar(cgf, 'eastoutside');


saveas(cgf, 'FVA_similarity_heatmap_H_.png');
saveas(cgf, 'FVA_similarity_heatmap_H_.pdf');