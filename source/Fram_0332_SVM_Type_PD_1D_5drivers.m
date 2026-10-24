% Plot 1D partial dependences for the trained SVM type model.
%
% calcPlot = true  -> calculate partial dependencies using the trained model
%                    and save results
% calcPlot = false -> load previously saved partial dependences and only plot

Fram_0310_SVM_Type_Init_5drivers;
calcPlot = true;

%% Parameters

nGrid = 21; %51
B = 50;     %50

% File names depend on LM_Model, which includes PredictorSetID from init.
modelFileName = fullfile( ...
    HGP.Data.Path2ResultsData, ...
    ['trainedLocationModel_' LM_Model '.mat']);

pdFileName = fullfile( ...
    HGP.Data.Path2ResultsData, ...
    ['PDP_1D_SVM_Type_' LM_Model '.mat']);

%% Load trained model

tmp = load(modelFileName, 'Mdl1', 'CVStats');
Mdl1 = tmp.Mdl1;
if isfield(tmp, 'CVStats')
    CVStats = tmp.CVStats;
else
    CVStats = struct();
end
clearvars tmp;

classNames = string(Mdl1.ClassNames);
K = numel(classNames);
ClassifID = CVStats.ClassifID;
bestParams = CVStats.BestParams;

predNames = string(Mdl1.PredictorNames);
predNamesPlot = renamePredictors(predNames);
Xpd = X_train(:, cellstr(predNames));

fun = @(X) ml_returnPosterior(Mdl1, X);

switch ClassifID
    case 'ECOC_SVM_Auto_Cost'
        classes = categories(Y_train);
        Kclass = numel(classes);
        GroupCounts = grpstats(table(Y_train, 'VariableNames', {'Y'}), 'Y', 'numel');
        GroupCounts.Weight = 1./GroupCounts.GroupCount;

        Cost = ones(Kclass) - eye(Kclass);
        for k = 1:Kclass
            Cost(k, :) = Cost(k, :) * GroupCounts.Weight(k);
        end

        Cost_Options = {'Cost', Cost};
        ml_train_func = @(X, Y, params, varargin) ml_train_ECOC(X, Y, params, Cost_Options{:}, varargin{:});

    otherwise
        ml_train_func = @ml_train_ECOC;
end

%% Predictor order

if isfield(CVStats, 'PermutationImportance')
    Importance = CVStats.PermutationImportance;
else
    Importance = permutationImportance(Mdl1);
end

[~, impOrder] = sort(Importance.ImportanceMean, 'descend');
importancePredNames = string(Importance.Predictor(impOrder));

[tf, predOrder] = ismember(importancePredNames, predNames);
predOrder = predOrder(tf);
predNames = predNames(predOrder);
predNamesPlot = predNamesPlot(predOrder);

nPred = numel(predNames);

%% Calculate or load partial dependences

if calcPlot
    pdMeanAll = cell(nPred, 1);
    pdLowAll = cell(nPred, 1);
    pdHighAll = cell(nPred, 1);
    xsAll = cell(nPred, 1);
    N = height(Xpd);

    for j = 1:nPred
        j
        predj = char(predNames(j));

        xj = Xpd.(predj);
        xs = linspace(min(xj), max(xj), nGrid);

        pdFull = partialDependence( ...
            fun, ...
            predj, ...
            Xpd, ...
            "UseParallel", true, ...
            "QueryPoints", xs', ...
            "OutputColumns", 1:K);

        pdBoot = NaN(B, K, nGrid);

        parfor b = 1:B
            b
            idxBoot = randi(N, N, 1);
            Xb = Xpd(idxBoot, :);
            Yb = Y_train(idxBoot);

            switch ClassifID
                case {'ECOC_SVM', 'ECOC_SVM_Auto_Cost'}
                    Mdl_b = ml_train_func(Xb, Yb, bestParams, 'FitPosterior', true);
                otherwise
                    Mdl_b = ml_train_func(Xb, Yb, bestParams);
            end

            funB = @(X) ml_returnPosterior(Mdl_b, X);
            pd_b = partialDependence( ...
                funB, ...
                predj, ...
                Xb, ...
                "UseParallel", false, ...
                "QueryPoints", xs', ...
                "OutputColumns", 1:K);

            pdBoot(b, :, :) = pd_b;
        end

        pdLow = squeeze(prctile(pdBoot, 2.5, 1));
        pdHigh = squeeze(prctile(pdBoot, 97.5, 1));

        pdMeanAll{j} = pdFull;
        pdLowAll{j} = pdLow;
        pdHighAll{j} = pdHigh;
        xsAll{j} = xs;
    end

    PD = struct();
    PD.pdMeanAll = pdMeanAll;
    PD.pdLowAll = pdLowAll;
    PD.pdHighAll = pdHighAll;
    PD.xsAll = xsAll;
    PD.predNames = predNames;
    PD.predNamesPlot = predNamesPlot;
    PD.classNames = classNames;
    PD.nGrid = nGrid;
    PD.B = B;
    PD.K = K;
    PD.LM_Model = LM_Model;
    PD.PredictorSetID = PredictorSetID;

    save(pdFileName, 'PD');
else
    tmp = load(pdFileName, 'PD');
    PD = tmp.PD;
    clearvars tmp;

    pdMeanAll = PD.pdMeanAll;
    pdLowAll = PD.pdLowAll;
    pdHighAll = PD.pdHighAll;
    xsAll = PD.xsAll;
    predNames = PD.predNames;
    predNamesPlot = PD.predNamesPlot;
    classNames = PD.classNames;
    nGrid = PD.nGrid;
    B = PD.B;
    K = PD.K;
    nPred = numel(predNames);
end

%% Plot

[fg, MaxFigure] = addFindFigure( ...
    'SVM type 1D partial dependences', ...
    MaxFigure, ...
    [10, 20, 1145, 656], ...
    'normal', ...
    'on');

tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
clrs = f_Clrs_fresh(K);

for j = 1:nPred
    pd = pdMeanAll{j};
    pdLow = pdLowAll{j};
    pdHigh = pdHighAll{j};
    xs = xsAll{j};

    pdRange = max(pd, [], 2) - min(pd, [], 2);
    lineWidths = 1 + 4 * pdRange;

    ax = nexttile;
    hold(ax, 'on');
    colororder(ax, clrs);

    for k = 1:K
        xPatch = [xs, fliplr(xs)];
        yPatch = [pdLow(k, :), fliplr(pdHigh(k, :))];

        p = patch('XData', xPatch, 'YData', yPatch, ...
            'FaceColor', clrs(k, :), ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none', ...
            'Parent', ax);
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';

        plot(ax, xs, pd(k, :), ...
            'LineWidth', lineWidths(k), ...
            'Color', clrs(k, :));
    end

    xlabel(ax, predNamesPlot(j), 'Interpreter', 'none');
    ylabel(ax, 'Partial dependence');
   % title(ax, sprintf('%s', predNamesPlot(j)), 'Interpreter', 'none');
    ylim([0, 0.8])
    hold(ax, 'off');
    drawnow;
end

lg = legend(classNames, 'NumColumns', min(3, K), 'Location', 'best');
lg.Title.String = 'Community types';
drawnow
AddLetters2Plots(fg, {'c', 'd', 'e', 'f', 'g', 'h'}, "HShift", 0)
