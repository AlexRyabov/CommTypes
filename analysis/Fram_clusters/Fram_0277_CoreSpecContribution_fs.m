% Fram_0277_CoreSpecContribution_fs.m
%
% Global ridge model with bootstrap feature stability selection.
% Output is limited to the selected-species contribution plot and test
% predicted-vs-observed diagnostics.

clearvars
proj_Initialize;

%% Load data

[Data, DI] = dt_loadData(HGP, []); %#ok<ASGLU> 

tStPredictors = Data.tStPredictors;
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));

coreDataFile = fullfile(HGP.Data.Path2ResultsData, 'CommTypeCoreSpeciesData.mat');

if ~isfile(coreDataFile)
    error(['Core species data file was not found: %s\n' ...
        'Run Fram_0275_CommDiffandOverlaps first.'], coreDataFile);
end

S = load(coreDataFile, 'coreSpeciesData');
coreSpeciesData = S.coreSpeciesData;

if ~ismember('ClustID', tDMC.Properties.VariableNames)
    error('Column ClustID was not found in tDMC.');
end

if ~ismember('chlorophyll_a_0', tStPredictors.Properties.VariableNames)
    error('Column chlorophyll_a_0 was not found in tStPredictors.');
end

nSamples = size(coreSpeciesData.Abundance, 1);

if height(tDMC) ~= nSamples
    error('Number of rows in tDMC does not match the core species abundance matrix.');
end

if height(tStPredictors) ~= nSamples
    error('Number of rows in tStPredictors does not match the core species abundance matrix.');
end

%% Parameters and metadata

Fram_0277_CoreSpecContribution_init;

nBootstrap = 100;
bootstrapQuantileProbs = [0.1 0.90];
minSelectedPredictors = 5;
fallbackTopN = 10;
finalModelMethod = "ridge"; % "ridge" or "fitlm"

speciesNames = string(coreSpeciesData.Species(:));
speciesPhyla = string(coreSpeciesData.Phylum(:));
speciesRoles = string(coreSpeciesData.Role(:));

if isfield(coreSpeciesData, 'PlotName')
    speciesPlotNames = string(coreSpeciesData.PlotName(:));
else
    speciesPlotNames = speciesNames;
end

predictorNames = matlab.lang.makeValidName(speciesNames);
predictorNames = matlab.lang.makeUniqueStrings(predictorNames);

Xraw = double(coreSpeciesData.Abundance);
coreSpeciesRowSum = sum(Xraw, 2, 'omitnan');
idxZeroCoreSpeciesRows = coreSpeciesRowSum <= 0 | isnan(coreSpeciesRowSum);

if any(idxZeroCoreSpeciesRows)
    warning('%d samples have zero or invalid total abundance across core species and will be excluded.', ...
        nnz(idxZeroCoreSpeciesRows));
end

Xraw(~idxZeroCoreSpeciesRows, :) = Xraw(~idxZeroCoreSpeciesRows, :) ./ coreSpeciesRowSum(~idxZeroCoreSpeciesRows);
Xraw(idxZeroCoreSpeciesRows, :) = NaN;

YoriginalAll = tStPredictors.chlorophyll_a_0;
clustIDAll = tDMC.ClustID(:);
communityTypes = sort(unique(clustIDAll(:)))';

switch modelAbundanceTransform
    case "clr"
        [XmodelAll, abundanceInfo] = local_clrTransformAbundance(Xraw);
    case "none"
        XmodelAll = Xraw;
        abundanceInfo = local_abundanceInfo("none", 0, NaN);
    case "power025"
        XmodelAll = local_powerTransformAbundance(Xraw, 0.25);
        abundanceInfo = local_abundanceInfo("power025", 0, NaN);
    case "power05"
        XmodelAll = local_powerTransformAbundance(Xraw, 0.5);
        abundanceInfo = local_abundanceInfo("power05", 0, NaN);
    otherwise
        error('Unknown modelAbundanceTransform: %s.', modelAbundanceTransform);
end

abundanceInfo.InputScale = "core_species_normalized_proportions";
abundanceInfo.CoreSpeciesRowNormalized = true;
abundanceInfo.NZeroCoreSpeciesRows = nnz(idxZeroCoreSpeciesRows);

validRows = ~isnan(YoriginalAll) & YoriginalAll > 0 & all(~isnan(XmodelAll), 2);

X = XmodelAll(validRows, :);
Yoriginal = YoriginalAll(validRows);
Y = chlTransform(Yoriginal);
clustID = clustIDAll(validRows);
cruiseID = tStPredictors.Cruise(validRows);
monthNum = tStPredictors.MonthNum(validRows);

%% Hold-out split

idxTest = false(size(Y));

for iType = 1:numel(communityTypes)
    idxType = find(clustID == communityTypes(iType));
    nTestType = max(1, round(holdoutFraction * numel(idxType)));
    idxType = idxType(randperm(numel(idxType)));
    idxTest(idxType(1:nTestType)) = true;
end

idxTrain = ~idxTest;

idxVariableSpecies = std(X(idxTrain, :), 0, 1, 'omitnan') > 0;
X = X(:, idxVariableSpecies);

modelPredictorNames = predictorNames(idxVariableSpecies);
modelSpeciesNames = speciesNames(idxVariableSpecies);
modelSpeciesPhyla = speciesPhyla(idxVariableSpecies);
modelSpeciesRoles = speciesRoles(idxVariableSpecies);
modelSpeciesPlotNames = speciesPlotNames(idxVariableSpecies);

nPredictors = size(X, 2);

%% First ridge pass: choose lambda once

[XzTrainAll, xMuAll, xSigmaAll] = local_standardizeTrain(X, idxTrain);
sampleWeights = local_clusterBalancedWeights(clustID, idxTrain, communityTypes);
nModelFolds = min(nFolds, sum(idxTrain));
cv_inner = [];

if blockedCV
    nModelFolds = min(nModelFolds, numel(unique(cruiseID(idxTrain))));
    [~, ~, ~, ~, cv_inner] = ml_blockedCV(XzTrainAll(idxTrain, :), Y(idxTrain), ...
        cruiseID(idxTrain), monthNum(idxTrain), nModelFolds);
    nModelFolds = cv_inner.NumTestSets;

end

[bestLambda, cvMSE] = local_selectRidgeLambda( ...
    XzTrainAll(idxTrain, :), Y(idxTrain), lambdaGrid, lambdaSelectionTolerance, ...
    cv_inner, nModelFolds, sampleWeights(idxTrain));

%% Bootstrap feature stability selection

betaBoot = nan(nPredictors, nBootstrap);
idxTrainByCluster = cell(numel(communityTypes), 1);
clusterTrainCounts = zeros(numel(communityTypes), 1);

for iType = 1:numel(communityTypes)
    idxTrainByCluster{iType} = find(clustID == communityTypes(iType) & idxTrain);
    clusterTrainCounts(iType) = numel(idxTrainByCluster{iType});
end

idxNonEmptyClusters = clusterTrainCounts > 0;
idxBootstrapClusters = find(idxNonEmptyClusters);
bootstrapSamplesPerCluster = max(clusterTrainCounts(idxNonEmptyClusters));

parfor iBoot = 1:nBootstrap
    idxBoot = zeros(bootstrapSamplesPerCluster * numel(idxBootstrapClusters), 1);

    for ii = 1:numel(idxBootstrapClusters)
        iType = idxBootstrapClusters(ii);
        idxTypeTrain = idxTrainByCluster{iType};
        idxBootType = idxTypeTrain(randi(numel(idxTypeTrain), bootstrapSamplesPerCluster, 1));
        idxBoot((ii-1)*bootstrapSamplesPerCluster + (1:bootstrapSamplesPerCluster)) = idxBootType(:);
    end

    xMuBoot = mean(X(idxBoot, :), 1, 'omitnan');
    xSigmaBoot = std(X(idxBoot, :), 0, 1, 'omitnan');
    xSigmaBoot(xSigmaBoot == 0 | isnan(xSigmaBoot)) = 1;
    XzBoot = (X(idxBoot, :) - xMuBoot) ./ xSigmaBoot;

    mdlBoot = fitrlinear(XzBoot, Y(idxBoot), ...
        'Learner', 'leastsquares', ...
        'Regularization', 'ridge', ...
        'Lambda', bestLambda);

    betaBoot(:, iBoot) = mdlBoot.Beta(:);
end

betaQ = quantile(betaBoot, bootstrapQuantileProbs, 2);
betaMedian = median(betaBoot, 2, 'omitnan');
idxSelected = betaQ(:, 1) > 0 | betaQ(:, 2) < 0;

if nnz(idxSelected) < minSelectedPredictors
    [~, idxFallback] = sort(abs(betaMedian), 'descend', 'MissingPlacement', 'last');
    idxSelected(idxFallback(1:min(fallbackTopN, numel(idxFallback)))) = true;
    warning('Bootstrap selected fewer than %d predictors. Falling back to top %d by |median beta|.', ...
        minSelectedPredictors, min(fallbackTopN, numel(idxFallback)));
end

fprintf('Bootstrap selected %d of %d predictors.\n', nnz(idxSelected), nPredictors);

%% Refit selected predictors

Xsel = X(:, idxSelected);
[XzSel, xMuSel, xSigmaSel] = local_standardizeTrain(Xsel, idxTrain);

bestLambdaSelected = NaN;
cvMSESelected = [];

switch finalModelMethod
    case "ridge"
        [bestLambdaSelected, cvMSESelected] = local_selectRidgeLambda( ...
            XzSel(idxTrain, :), Y(idxTrain), lambdaGrid, lambdaSelectionTolerance, ...
            cv_inner, nModelFolds, sampleWeights(idxTrain));

        linearMdl = fitrlinear(XzSel(idxTrain, :), Y(idxTrain), ...
            'Learner', 'leastsquares', ...
            'Regularization', 'ridge', ...
            'Lambda', bestLambdaSelected, ...
            'Weights', sampleWeights(idxTrain));

        YhatTransformed = predict(linearMdl, XzSel);
        betaStandardized = linearMdl.Beta(:);

    case "fitlm"
        predictorNamesSelected = matlab.lang.makeValidName(modelSpeciesNames(idxSelected));
        predictorNamesSelected = matlab.lang.makeUniqueStrings(predictorNamesSelected);
        tFit = array2table(XzSel(idxTrain, :), 'VariableNames', cellstr(predictorNamesSelected));
        tFit.ChlorophyllTransformed = Y(idxTrain);

        linearMdl = fitlm(tFit, 'linear', ...
            'ResponseVar', 'ChlorophyllTransformed', ...
            'PredictorVars', cellstr(predictorNamesSelected), ...
            'Weights', sampleWeights(idxTrain));

        tPredict = array2table(XzSel, 'VariableNames', cellstr(predictorNamesSelected));
        YhatTransformed = predict(linearMdl, tPredict);
        coefTbl = linearMdl.Coefficients;
        coefNames = string(coefTbl.Properties.RowNames);
        betaStandardized = zeros(size(XzSel, 2), 1);
        [isCoef, idxCoef] = ismember(string(predictorNamesSelected), coefNames);
        betaStandardized(isCoef) = coefTbl.Estimate(idxCoef(isCoef));

    otherwise
        error('Unknown finalModelMethod: %s.', finalModelMethod);
end

YhatOriginal = chlBackTransform(YhatTransformed);
betaUnstandardized = betaStandardized ./ xSigmaSel(:);

modelSpeciesNames = modelSpeciesNames(idxSelected);
modelSpeciesPhyla = modelSpeciesPhyla(idxSelected);
modelSpeciesRoles = modelSpeciesRoles(idxSelected);
modelSpeciesPlotNames = modelSpeciesPlotNames(idxSelected);

%% Plot selected-predictor ridge lambda CV loss

if finalModelMethod == "ridge" && ~isempty(cvMSESelected)
    [fgLoss, MaxFigure] = addFindFigure('FS ridge selected predictors lambda CV loss', ...
        MaxFigure, [40, 80, 720, 470], 'normal', 'on');
    clf(fgLoss);
    fgLoss.Color = 'w';

    axLoss = axes(fgLoss);
    semilogx(axLoss, lambdaGrid, cvMSESelected, '-o', ...
        'Color', [0.20 0.38 0.80], ...
        'MarkerFaceColor', [0.20 0.38 0.80], ...
        'MarkerSize', 4, ...
        'LineWidth', 1.2);
    hold(axLoss, 'on');
    xline(axLoss, bestLambdaSelected, 'k--', ...
        sprintf('\\lambda = %.3g', bestLambdaSelected), ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom');
    grid(axLoss, 'on');
    xlabel(axLoss, '\lambda');
    ylabel(axLoss, 'CV MSE');
    title(axLoss, sprintf('Selected predictors ridge CV loss, n = %d', nnz(idxSelected)));
end

%% Contribution summaries by community type

nSelectedPredictors = nnz(idxSelected);
medianContribution = nan(nSelectedPredictors, numel(communityTypes));
contributionQLow = nan(nSelectedPredictors, numel(communityTypes));
contributionQHigh = nan(nSelectedPredictors, numel(communityTypes));

for iType = 1:numel(communityTypes)
    idxType = clustID == communityTypes(iType);
    sampleContribType = Xsel(idxType, :) .* betaUnstandardized(:)';

    for iSpecies = 1:nSelectedPredictors
        vals = sampleContribType(:, iSpecies);
        vals = vals(~isnan(vals));

        if isempty(vals)
            continue
        end

        Q = quantile(vals, contributionQuantileProbs);
        contributionQLow(iSpecies, iType) = Q(1);
        medianContribution(iSpecies, iType) = Q(2);
        contributionQHigh(iSpecies, iType) = Q(3);
    end
end

%% Plot selected core species contributions

[fgB, MaxFigure] = addFindFigure('FS ridge selected species contributions', ...
    MaxFigure, [24 332 1461 455], 'normal', 'on');

if finalModelMethod == "ridge"
    ttl = sprintf('Bootstrap-selected ridge species contribution (n = %d, lambda = %.3g)', ...
        nSelectedPredictors, bestLambdaSelected);
else
    ttl = sprintf('Bootstrap-selected fitlm species contribution (n = %d)', nSelectedPredictors);
end

barOpts = struct();
barOpts.Figure = fgB;
barOpts.ClearFigure = true;
barOpts.XLabel = 'Community type';
barOpts.YLabel = char(coreSpeciesData.GroupName);
barOpts.Title = ttl;
barOpts.ColorbarLabel = 'median unstandardized predictor \times beta';
barOpts.RoleColors = roleClrs;
barOpts.RoleShortNames = roleShortNames;
barOpts.RoleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};
barOpts.FaceAlpha = 0.88;
barOpts.XGridLineColor = [0.55 0.55 0.55];

typeLabels = strings(size(communityTypes));

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    idxTypeTest = clustID == typeID & idxTest;
    idxValid = idxTypeTest & Yoriginal > 0 & YhatOriginal > 0;
    r2LogType = hf_stat_R2(log10(Yoriginal(idxValid)), log10(YhatOriginal(idxValid)));
    rCorr = corr(log10(Yoriginal(idxValid)), log10(YhatOriginal(idxValid)), 'Rows', 'complete');
    typeLabels(iType) = sprintf('Type %d \\newlineR^2_{log}=%.2f \\newline r_{log}=%.2f', typeID, r2LogType, rCorr);
end

[fgB, axB, cb, maxAbsPlotContribution] = plotGroupedIntervalBarsHor( ...
    medianContribution, contributionQLow, contributionQHigh, ...
    typeLabels, modelSpeciesPlotNames, ...
    modelSpeciesPhyla, modelSpeciesRoles, barOpts);

%% Test predicted/observed diagnostics

tPred = table();
tPred.ClustID = clustID;
tPred.ChlorophyllObserved = Yoriginal;
tPred.ChlorophyllPredicted = YhatOriginal;
tPred.DataSet = strings(size(Y));
tPred.DataSet(idxTrain) = "train";
tPred.DataSet(idxTest) = "test";

[fgPred, MaxFigure] = addFindFigure('FS ridge selected species observed predicted test', ...
    MaxFigure, [10, 20, 1100 650], 'normal', 'on');
clf(fgPred);
fgPred.Color = 'w';

tlPred = tiledlayout(fgPred, 'flow', 'Padding', 'compact', 'TileSpacing', 'compact');
idxSet = tPred.DataSet == "test";
predAxisLim = local_positivePlotLimits(tPred.ChlorophyllObserved(idxSet), tPred.ChlorophyllPredicted(idxSet));

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    axPred = nexttile(tlPred);
    idxType = tPred.ClustID == typeID & idxSet;
    idxPlot = idxType & tPred.ChlorophyllObserved > 0 & tPred.ChlorophyllPredicted > 0;

    scatter(axPred, tPred.ChlorophyllPredicted(idxPlot), ...
        tPred.ChlorophyllObserved(idxPlot), 22, ...
        'filled', ...
        'MarkerFaceColor', [0.20 0.38 0.80], ...
        'MarkerFaceAlpha', 0.65, ...
        'MarkerEdgeColor', [0.25 0.25 0.25]);
    hold(axPred, 'on');

    if all(isfinite(predAxisLim)) && predAxisLim(2) > predAxisLim(1)
        plot(axPred, predAxisLim, predAxisLim, 'k-', 'LineWidth', 0.75);
    end

    yObsO = tPred.ChlorophyllObserved(idxType);
    yHatO = tPred.ChlorophyllPredicted(idxType);
    r2Log = hf_stat_R2(log10(yObsO), log10(yHatO));
    rLog = local_corrValue(log10(yObsO), log10(yHatO));

    set(axPred, 'XScale', 'log', 'YScale', 'log');
    xlim(axPred, predAxisLim);
    ylim(axPred, predAxisLim);
    axis(axPred, 'square');
    grid(axPred, 'on');
    xlabel(axPred, 'Predicted chlorophyll a');
    ylabel(axPred, 'Observed chlorophyll a');
    title(axPred, sprintf('test Type %d, R^2_{log}=%.2f, r_{log}=%.2f', ...
        typeID, r2Log, rLog));
end

disp('Finished.');

function info = local_abundanceInfo(transformName, pseudocount, maxAbsRowSum)

    info = struct();
    info.Transform = transformName;
    info.InputScale = "proportions";
    info.Pseudocount = pseudocount;
    info.MaxAbsRowSum = maxAbsRowSum;
end

function X = local_powerTransformAbundance(X, alpha)

    X = double(X);

    if any(X(:) < 0)
        error('Power transform requires non-negative abundance values.');
    end

    X = X .^ alpha;
end

function [Xclr, info] = local_clrTransformAbundance(X)

    X = double(X);
    pseudocount = 1e-6;

    if any(X(:) < 0)
        error('CLR transform requires non-negative abundance values.');
    end

    X = X + pseudocount;
    rowSums = sum(X, 2, 'omitnan');

    if any(rowSums <= 0 | isnan(rowSums))
        error('CLR transform found rows with non-positive abundance sums.');
    end

    X = X ./ rowSums;
    logX = log10(X);
    Xclr = logX - mean(logX, 2, 'omitnan');

    info = local_abundanceInfo("clr", pseudocount, max(abs(sum(Xclr, 2, 'omitnan')), [], 'omitnan'));
end

function [bestLambda, cvMSE] = local_selectRidgeLambda(Xtrain, Ytrain, lambdaGrid, lambdaSelectionTolerance, cv_inner, nModelFolds, sampleWeights)

    if ~isempty(cv_inner)
        cvMdl = fitrlinear(Xtrain, Ytrain, ...
            'Learner', 'leastsquares', ...
            'Regularization', 'ridge', ...
            'CVPartition', cv_inner, ...
            'Lambda', lambdaGrid, ...
            'Weights', sampleWeights);
    else
        cvMdl = fitrlinear(Xtrain, Ytrain, ...
            'Learner', 'leastsquares', ...
            'Regularization', 'ridge', ...
            'KFold', nModelFolds, ...
            'Lambda', lambdaGrid, ...
            'Weights', sampleWeights);
    end

    cvMSE = kfoldLoss(cvMdl);
    minCVMSE = min(cvMSE);
    idxAllowedLambda = cvMSE <= minCVMSE * (1 + lambdaSelectionTolerance);
    idxBestLambda = find(idxAllowedLambda, 1, 'last');
    bestLambda = lambdaGrid(idxBestLambda);
end

function [Xz, xMu, xSigma] = local_standardizeTrain(X, idxTrain)

    xMu = mean(X(idxTrain, :), 1, 'omitnan');
    xSigma = std(X(idxTrain, :), 0, 1, 'omitnan');
    xSigma(xSigma == 0 | isnan(xSigma)) = 1;
    Xz = (X - xMu) ./ xSigma;
end

function sampleWeights = local_clusterBalancedWeights(clustID, idxTrain, communityTypes)

    sampleWeights = ones(size(clustID));
    sampleWeights(:) = 0;

    for iType = 1:numel(communityTypes)
        idxTypeTrain = clustID == communityTypes(iType) & idxTrain;

        if any(idxTypeTrain)
            sampleWeights(idxTypeTrain) = 1 / sum(idxTypeTrain);
        end
    end

    sampleWeights(idxTrain) = sampleWeights(idxTrain) ./ mean(sampleWeights(idxTrain), 'omitnan');
    sampleWeights(~idxTrain) = 1;
end

function labels = local_coloredYTickLabels(names, roles, roleClrs)

    names = string(names);
    roles = string(roles);
    labels = strings(size(names));

    for ii = 1:numel(names)
        clr = local_roleColor(roles(ii), roleClrs);
        labels(ii) = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%s', ...
            clr(1), clr(2), clr(3), local_texEscape(names(ii)));
    end
end

function out = local_texEscape(in)

    out = char(in);
    out = strrep(out, '\', '\\');
    out = strrep(out, '_', '\_');
    out = strrep(out, '^', '\^');
    out = strrep(out, '{', '\{');
    out = strrep(out, '}', '\}');
end

function clr = local_roleColor(roleName, roleClrs)

    roleName = char(roleName);

    if isKey(roleClrs, roleName)
        clr = roleClrs(roleName);
    else
        clr = roleClrs('Unknown');
    end
end

function local_addRoleLegend(ax, roleClrs, roleShortNames, fontSize)

    roles = ["Autotroph", "Mixotroph", "Zooplankton", "Parasite", "Recycler"];
    x0 = -0.13;
    dx = 0.075;
    y0 = 1.055;

    for ii = 1:numel(roles)
        clr = local_roleColor(roles(ii), roleClrs);
        text(ax, x0 + (ii-1)*dx, y0, roleShortNames(char(roles(ii))), ...
            'Units', 'normalized', ...
            'Color', clr, ...
            'FontWeight', 'bold', ...
            'FontSize', fontSize, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'Clipping', 'off');
    end
end

function cmap = local_blueWhiteRed(n)

    if nargin < 1
        n = 256;
    end

    n1 = floor(n/2);
    n2 = n - n1;
    blue = [0.18 0.38 0.82];
    white = [1 1 1];
    red = [0.84 0.18 0.12];
    cmap1 = [linspace(blue(1), white(1), n1)', ...
        linspace(blue(2), white(2), n1)', ...
        linspace(blue(3), white(3), n1)'];
    cmap2 = [linspace(white(1), red(1), n2)', ...
        linspace(white(2), red(2), n2)', ...
        linspace(white(3), red(3), n2)'];
    cmap = [cmap1; cmap2];
end

function clr = local_valueColor(value, maxAbsValue, cmap)

    if isnan(value) || maxAbsValue <= 0 || isnan(maxAbsValue)
        clr = [0.85 0.85 0.85];
        return
    end

    x = max(-1, min(1, value / maxAbsValue));
    idx = 1 + round((x + 1) / 2 * (size(cmap, 1) - 1));
    clr = cmap(idx, :);
end

function r = local_corrValue(y, yHat)

    idx = ~isnan(y) & ~isnan(yHat) & isfinite(y) & isfinite(yHat);
    y = y(idx);
    yHat = yHat(idx);

    if numel(y) < 2
        r = NaN;
        return
    end

    r = corr(y(:), yHat(:), 'Rows', 'complete');
end

function lim = local_positivePlotLimits(varargin)

    vals = [];

    for ii = 1:nargin
        x = varargin{ii};
        vals = [vals; x(:)]; %#ok<AGROW>
    end

    vals = vals(vals > 0 & isfinite(vals));

    if isempty(vals)
        lim = [1e-3 1];
        return
    end

    lim = [min(vals) max(vals)];

    if lim(1) == lim(2)
        lim = lim .* [0.8 1.2];
    end

    lim = 10 .^ (log10(lim) + [-0.08 0.08]);
end
