% Fram_0296_CoreContribFinalClean.m
%
% Clean selected-taxa workflow:
% load all core indicator taxa, row-normalize abundances, transform predictors,
% use stepwise selection on the full normalized predictor set, then rebuild
% predictors from raw abundances for selected taxa only and fit the final
% linear chlorophyll model.

clearvars
proj_Initialize;

%% Load data

[Data, ~] = dt_loadData(HGP, []);

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

if ~ismember('MonthNum', tStPredictors.Properties.VariableNames)
    error('Column MonthNum was not found in tStPredictors.');
end

if ~ismember('Cruise', tStPredictors.Properties.VariableNames)
    error('Column Cruise was not found in tStPredictors.');
end

nSamples = size(coreSpeciesData.Abundance, 1);

if height(tDMC) ~= nSamples
    error('Number of rows in tDMC does not match the core taxa abundance matrix.');
end

if height(tStPredictors) ~= nSamples
    error('Number of rows in tStPredictors does not match the core taxa abundance matrix.');
end

%% Initialize settings and metadata

Fram_0290_CoreSpecContribution_init;

selectionModelMethod = "stepwise";
finalSelectedModelMethod = "fitlm";
speciesNames = string(coreSpeciesData.Species(:));
speciesPhyla = string(coreSpeciesData.Phylum(:));
speciesRoles = string(coreSpeciesData.Role(:));

if isfield(coreSpeciesData, 'PlotName')
    speciesPlotNames = string(coreSpeciesData.PlotName(:));
else
    speciesPlotNames = speciesNames;
end

%% Transform and normalize core indicator taxa abundances

Xraw = double(coreSpeciesData.Abundance);
coreTaxaRowSum = sum(Xraw, 2, 'omitnan');
idxZeroCoreTaxaRows = coreTaxaRowSum <= 0 | isnan(coreTaxaRowSum);

if any(idxZeroCoreTaxaRows)
    warning('%d samples have zero or invalid total abundance across core taxa and will be excluded.', ...
        nnz(idxZeroCoreTaxaRows));
end

Xraw(idxZeroCoreTaxaRows, :) = NaN;

XmodelAll = local_normalizePositiveRows(Xraw);
XmodelAll = local_transformAbundance(XmodelAll, modelAbundanceTransform, modelAbundanceLogPseudocount);

YoriginalAll = tStPredictors.chlorophyll_a_0;
YtransformedAll = chlTransform(YoriginalAll);
clustIDAll = tDMC.ClustID(:);
monthNumAll = tStPredictors.MonthNum(:);
cruiseIDAll = tStPredictors.Cruise(:);

idxExcludedLowChlorophyllTypes = false(size(YoriginalAll));

if excludeLowChlorophyllTypes
    idxExcludedLowChlorophyllTypes = ismember(clustIDAll, lowChlorophyllExcludedTypes) & ...
        YoriginalAll < lowChlorophyllExclusionThreshold;

    if any(idxExcludedLowChlorophyllTypes)
        fprintf('Excluded %d samples with ClustID in [%s] and chlorophyll_a_0 < %.4g.\n', ...
            nnz(idxExcludedLowChlorophyllTypes), ...
            strjoin(string(lowChlorophyllExcludedTypes), ', '), ...
            lowChlorophyllExclusionThreshold);
    end
end

validRows = ~isnan(YoriginalAll) & YoriginalAll > 0 & ...
    ~isnan(YtransformedAll) & all(~isnan(XmodelAll), 2) & ...
    ~isnan(clustIDAll) & ~isnan(monthNumAll) & ...
    ~idxExcludedLowChlorophyllTypes;

if ~any(validRows)
    error('No valid rows remain after applying data filters.');
end

X = XmodelAll(validRows, :);
Y = YtransformedAll(validRows);
Yoriginal = YoriginalAll(validRows);
clustID = clustIDAll(validRows);
monthNum = monthNumAll(validRows);
cruiseID = cruiseIDAll(validRows);
communityTypes = sort(unique(clustID(:)))';

%% Hold-out split

idxTest = false(size(Y));

for iType = 1:numel(communityTypes)
    idxType = find(clustID == communityTypes(iType));
    nTestType = max(1, round(holdoutFraction * numel(idxType)));
    idxType = idxType(randperm(numel(idxType)));
    idxTest(idxType(1:nTestType)) = true;
end

idxTrain = ~idxTest;

idxVariableTaxa = std(X(idxTrain, :), 0, 1, 'omitnan') > 0;
X = X(:, idxVariableTaxa);
modelSpeciesNames = speciesNames(idxVariableTaxa);
modelSpeciesPhyla = speciesPhyla(idxVariableTaxa);
modelSpeciesRoles = speciesRoles(idxVariableTaxa);
modelSpeciesPlotNames = speciesPlotNames(idxVariableTaxa);

%% Stepwise selection on all core taxa

[Xz, ~, ~] = local_standardizeTrain(X, idxTrain);
sampleWeights = local_clusterBalancedWeights(clustID, idxTrain, communityTypes, useClusterBalancedRidgeWeights);
modelPredictorNames = matlab.lang.makeValidName(modelSpeciesNames);
modelPredictorNames = matlab.lang.makeUniqueStrings(modelPredictorNames);

if selectionModelMethod ~= "stepwise"
    error('Unknown selectionModelMethod: %s.', selectionModelMethod);
end

tSelectionFit = array2table(Xz(idxTrain, :), 'VariableNames', cellstr(modelPredictorNames));
tSelectionFit.ChlorophyllTransformed = Y(idxTrain);

mdlSelection = stepwiselm(tSelectionFit, 'constant', ...
    'Upper', 'linear', ...
    'ResponseVar', 'ChlorophyllTransformed', ...
    'PredictorVars', cellstr(modelPredictorNames), ...
    'Criterion', char(stepwiseCriterion), ...
    'PEnter', stepwisePEnter, ...
    'PRemove', stepwisePRemove, ...
    'Verbose', 0);

selectionCoefNames = local_tableRowNames(mdlSelection.Coefficients);
idxSelectedModelTaxa = ismember(modelPredictorNames, selectionCoefNames);

if ~any(idxSelectedModelTaxa)
    error('Stepwise selection did not keep any taxa predictors.');
end

fprintf('Stepwise selection kept %d of %d variable core taxa predictors.\n', ...
    nnz(idxSelectedModelTaxa), numel(idxSelectedModelTaxa));

selectedTaxa = modelSpeciesNames(idxSelectedModelTaxa);
idxSelectedOriginal = ismember(speciesNames, selectedTaxa);

%% Rebuild selected-taxa predictors from raw abundances and fit final linear model

XselectedAll = local_normalizePositiveRows(Xraw(:, idxSelectedOriginal));
XselectedAll = local_transformAbundance(XselectedAll, ...
    modelAbundanceTransform, modelAbundanceLogPseudocount);
Xselected = XselectedAll(validRows, :);

selectedSpeciesNames = speciesNames(idxSelectedOriginal);
selectedSpeciesPhyla = speciesPhyla(idxSelectedOriginal);
selectedSpeciesRoles = speciesRoles(idxSelectedOriginal);
selectedSpeciesPlotNames = speciesPlotNames(idxSelectedOriginal);

idxVariableSelected = std(Xselected(idxTrain, :), 0, 1, 'omitnan') > 0;

if ~all(idxVariableSelected)
    warning('Final selected-taxa predictor matrix has %d constant predictors after selected-taxa renormalization.', ...
        nnz(~idxVariableSelected));
end

Xselected = Xselected(:, idxVariableSelected);
selectedSpeciesNames = selectedSpeciesNames(idxVariableSelected);
selectedSpeciesPhyla = selectedSpeciesPhyla(idxVariableSelected);
selectedSpeciesRoles = selectedSpeciesRoles(idxVariableSelected);
selectedSpeciesPlotNames = selectedSpeciesPlotNames(idxVariableSelected);

[XzSelected, ~, xSigmaSelected] = local_standardizeTrain(Xselected, idxTrain);
selectedPredictorNames = matlab.lang.makeValidName(selectedSpeciesNames);
selectedPredictorNames = matlab.lang.makeUniqueStrings(selectedPredictorNames);

tFit = array2table(XzSelected(idxTrain, :), 'VariableNames', cellstr(selectedPredictorNames));
tFit.ChlorophyllTransformed = Y(idxTrain);

mdlFinal = fitlm(tFit, 'linear', ...
    'ResponseVar', 'ChlorophyllTransformed', ...
    'PredictorVars', cellstr(selectedPredictorNames), ...
    'Weights', sampleWeights(idxTrain));

tPredict = array2table(XzSelected, 'VariableNames', cellstr(selectedPredictorNames));
YhatSelectedTransformed = predict(mdlFinal, tPredict);
YhatSelectedOriginal = chlBackTransform(YhatSelectedTransformed);

tFinalModelCoefficients = local_linearModelCoefficientTable(mdlFinal);
finalModelCoefficientFile = fullfile(HGP.Data.Path2ResultsData, ...
    sprintf('chl_predict_clean_%s_model.xlsx', finalSelectedModelMethod));
writetable(tFinalModelCoefficients, finalModelCoefficientFile);
fprintf('Saved clean final linear model coefficients to %s.\n', finalModelCoefficientFile);

coefTbl = mdlFinal.Coefficients;
coefNames = local_tableRowNames(coefTbl);
betaSelectedStandardized = zeros(size(XzSelected, 2), 1);
[isCoef, idxCoef] = ismember(string(selectedPredictorNames), coefNames);
betaSelectedStandardized(isCoef) = coefTbl.Estimate(idxCoef(isCoef));
betaSelectedUnstandardized = betaSelectedStandardized ./ xSigmaSelected(:);

fprintf('Final clean linear model fit: %d selected taxa predictors.\n', numel(selectedSpeciesNames));

%% Contribution by community type for final selected-taxa model

nSelectedTaxa = numel(selectedSpeciesNames);
nTypes = numel(communityTypes);
meanContributionByType = nan(nSelectedTaxa, nTypes);
contributionQLow = nan(nSelectedTaxa, nTypes);
contributionQHigh = nan(nSelectedTaxa, nTypes);

for iType = 1:nTypes
    idxType = clustID == communityTypes(iType);
    sampleContribType = Xselected(idxType, :) .* betaSelectedUnstandardized(:)';

    for iTaxon = 1:nSelectedTaxa
        vals = sampleContribType(:, iTaxon);
        vals = vals(~isnan(vals));

        if isempty(vals)
            continue
        end

        Q = quantile(vals, contributionQuantileProbs);
        contributionQLow(iTaxon, iType) = Q(1);
        meanContributionByType(iTaxon, iType) = mean(vals, 'omitnan');
        contributionQHigh(iTaxon, iType) = Q(3);
    end
end

[fgB, MaxFigure] = addFindFigure('Fram 0296 clean selected core taxa contributions by type', ...
    MaxFigure, [24 332 1131 455], 'normal', 'on');

barOpts = struct();
barOpts.Figure = fgB;
barOpts.ClearFigure = true;
barOpts.XLabel = 'Community type';
barOpts.YLabel = 'Indicator taxa';
barOpts.ColorbarLabel = 'mean transformed selected-taxa predictor \times beta';
barOpts.RoleColors = roleClrs;
barOpts.RoleShortNames = roleShortNames;
barOpts.RoleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};
barOpts.FaceAlpha = 0.88;
barOpts.XGridLineColor = [0.55 0.55 0.55];

typeLabels = strings(size(communityTypes));

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    idxTypeTest = clustID == typeID & idxTest;
    idxValid = idxTypeTest & isfinite(Y) & isfinite(YhatSelectedTransformed);
    r2TransfType = hf_stat_R2(Y(idxValid), YhatSelectedTransformed(idxValid));
    rTransfType = local_corrValue(Y(idxValid), YhatSelectedTransformed(idxValid));
    typeLabels(iType) = sprintf('Type %d \\newlineR^2_{transf}=%.2f \\newline r_{transf}=%.2f', ...
        typeID, r2TransfType, rTransfType);
end

plotGroupedIntervalBarsHor( ...
    meanContributionByType, contributionQLow, contributionQHigh, ...
    typeLabels, selectedSpeciesPlotNames, ...
    selectedSpeciesPhyla, selectedSpeciesRoles, barOpts);

%% Final selected-taxa predicted/observed diagnostics

tPredSelected = table();
tPredSelected.ClustID = clustID;
tPredSelected.ChlorophyllObserved = Yoriginal;
tPredSelected.ChlorophyllPredicted = YhatSelectedOriginal;
tPredSelected.ChlorophyllObservedTransformed = Y;
tPredSelected.ChlorophyllPredictedTransformed = YhatSelectedTransformed;
tPredSelected.DataSet = strings(size(Y));
tPredSelected.DataSet(idxTrain) = "train";
tPredSelected.DataSet(idxTest) = "test";

MaxFigure = local_plotPredictedObservedByType(tPredSelected, communityTypes, "test", MaxFigure);
MaxFigure = local_plotPredictedObservedByType(tPredSelected, communityTypes, "train", MaxFigure);

disp('Finished.');

function X = local_transformAbundance(X, transformName, logPseudocount)

    X = double(X);

    if any(X(:) < 0)
        error('Abundance transform requires non-negative abundance values.');
    end

    switch string(transformName)
        case "none"
            return
        case "log10"
            if logPseudocount <= 0 || ~isfinite(logPseudocount)
                error('modelAbundanceLogPseudocount must be a positive finite value.');
            end

            X = log10(X + logPseudocount);
        case "power025"
            X = X .^ 0.25;
        case "power05"
            X = sqrt(X);
        otherwise
            error('Unknown modelAbundanceTransform: %s.', transformName);
    end
end

function X = local_normalizePositiveRows(X)

    idxInvalid = any(isnan(X), 2);
    rowSums = sum(X, 2, 'omitnan');
    idxPositive = ~idxInvalid & rowSums > 0 & isfinite(rowSums);
    idxZero = ~idxInvalid & ~idxPositive;

    X(idxPositive, :) = X(idxPositive, :) ./ rowSums(idxPositive);
    X(idxZero, :) = 0;
    X(idxInvalid, :) = NaN;
end

function [Xz, xMu, xSigma] = local_standardizeTrain(X, idxTrain)

    xMu = mean(X(idxTrain, :), 1, 'omitnan');
    xSigma = std(X(idxTrain, :), 0, 1, 'omitnan');
    xSigma(xSigma == 0 | isnan(xSigma)) = 1;
    Xz = (X - xMu) ./ xSigma;
end

function tCoef = local_linearModelCoefficientTable(mdl)

    coefTbl = mdl.Coefficients;
    coefVarNames = local_tableVariableNames(coefTbl);
    predictor = local_tableRowNames(coefTbl);
    nCoef = size(coefTbl, 1);
    lowerCI = nan(nCoef, 1);
    upperCI = nan(nCoef, 1);

    if ismember("Lower", coefVarNames) && ismember("Upper", coefVarNames)
        lowerCI = coefTbl.Lower;
        upperCI = coefTbl.Upper;
    else
        ci = coefCI(mdl);

        if isnumeric(ci) && size(ci, 2) >= 2
            lowerCI = ci(:, 1);
            upperCI = ci(:, 2);
        end
    end

    tCoef = table();
    tCoef.Predictor = predictor;
    tCoef.Estimate = coefTbl.Estimate;
    tCoef.SE = coefTbl.SE;
    tCoef.pValue = coefTbl.pValue;
    tCoef.LowerCI = lowerCI;
    tCoef.UpperCI = upperCI;
end

function varNames = local_tableVariableNames(T)

    try
        varNames = string(T.Properties.VariableNames);
    catch
        varNames = string(T.Properties.VarNames);
    end
end

function rowNames = local_tableRowNames(T)

    try
        rowNames = string(T.Properties.RowNames);
    catch
        rowNames = string(T.Properties.ObsNames);
    end
end

function sampleWeights = local_clusterBalancedWeights(clustID, idxTrain, communityTypes, useClusterBalancedRidgeWeights)

    sampleWeights = ones(size(clustID));

    if ~useClusterBalancedRidgeWeights
        return
    end

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

function MaxFigure = local_plotPredictedObservedByType(tPredSelected, communityTypes, plotSet, MaxFigure)

    figureName = ['Fram 0296 clean selected taxa observed predicted ' char(plotSet)];
    [fgPred, MaxFigure] = addFindFigure(figureName, ...
        MaxFigure, [10, 20, 1100 650], 'normal', 'on');
    clf(fgPred);
    fgPred.Color = 'w';

    tlPred = tiledlayout(fgPred, 'flow', 'Padding', 'compact', 'TileSpacing', 'compact');
    idxSet = tPredSelected.DataSet == plotSet;
    predAxisLim = local_positivePlotLimits(tPredSelected.ChlorophyllObserved(idxSet), ...
        tPredSelected.ChlorophyllPredicted(idxSet));

    for iType = 1:numel(communityTypes)
        typeID = communityTypes(iType);
        axPred = nexttile(tlPred);
        idxType = tPredSelected.ClustID == typeID & idxSet;
        idxPlot = idxType & tPredSelected.ChlorophyllObserved > 0 & ...
            tPredSelected.ChlorophyllPredicted > 0;

        scatter(axPred, tPredSelected.ChlorophyllPredicted(idxPlot), ...
            tPredSelected.ChlorophyllObserved(idxPlot), 22, ...
            'filled', ...
            'MarkerFaceColor', [0.20 0.38 0.80], ...
            'MarkerFaceAlpha', 0.65, ...
            'MarkerEdgeColor', [0.25 0.25 0.25]);
        hold(axPred, 'on');

        if all(isfinite(predAxisLim)) && predAxisLim(2) > predAxisLim(1)
            plot(axPred, predAxisLim, predAxisLim, 'k-', 'LineWidth', 0.75);
        end

        yObsTransf = tPredSelected.ChlorophyllObservedTransformed(idxType);
        yHatTransf = tPredSelected.ChlorophyllPredictedTransformed(idxType);
        r2Transf = hf_stat_R2(yObsTransf, yHatTransf);
        rTransf = local_corrValue(yObsTransf, yHatTransf);

        set(axPred, 'XScale', 'log', 'YScale', 'log');
        xlim(axPred, predAxisLim);
        ylim(axPred, predAxisLim);
        axis(axPred, 'square');
        grid(axPred, 'on');
        xlabel(axPred, 'Predicted chlorophyll a');
        ylabel(axPred, 'Observed chlorophyll a');
        title(axPred, sprintf('%s Type %d, R^2_{transf}=%.2f, r_{transf}=%.2f', ...
            char(plotSet), typeID, r2Transf, rTransf));
    end
end

function lim = local_positivePlotLimits(varargin)

    vals = [];

    for iArg = 1:nargin
        x = varargin{iArg};
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
