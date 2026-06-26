% Fram_0277_CoreSpecContribution.m
%
% Fit one global stepwise linear model for chlorophyll using all samples.
% Then estimate the contribution of each core species in each community type as
% transformed abundance multiplied by the unstandardized global regression beta.

clearvars
proj_Initialize;

%% Load data

[Data, DI] = dt_loadData(HGP, []);

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

if ~isempty(DI) && isfield(DI, 'idxSamplCols')
    nSampleCols = numel(DI.idxSamplCols);
    if nSampleCols ~= nSamples
        warning('DI.idxSamplCols has %d sample columns, but core data has %d samples.', nSampleCols, nSamples);
    end
end

%% Parameters and metadata

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

Fram_0277_CoreSpecContribution_init;

%% Transform abundance once and fit one global model

switch modelAbundanceTransform
    case "clr"
        [XmodelAll, abundanceInfo] = local_clrTransformAbundance(Xraw);

        fprintf('CLR transform: input=%s, pseudocount=%g, max(abs(row sum))=%g\n', ...
            char(abundanceInfo.InputScale), abundanceInfo.Pseudocount, abundanceInfo.MaxAbsRowSum);

        if abundanceInfo.MaxAbsRowSum > 1e-8
            warning('CLR row sums are not close to zero: max(abs(row sum)) = %g.', abundanceInfo.MaxAbsRowSum);
        end

    case "none"
        XmodelAll = Xraw;
        abundanceInfo = struct();
        abundanceInfo.Transform = "none";
        abundanceInfo.InputScale = "proportions";
        abundanceInfo.Pseudocount = 0;
        abundanceInfo.MaxAbsRowSum = NaN;

    case "power025"
        XmodelAll = local_powerTransformAbundance(Xraw, 0.25);
        abundanceInfo = struct();
        abundanceInfo.Transform = "power025";
        abundanceInfo.InputScale = "proportions";
        abundanceInfo.Pseudocount = 0;
        abundanceInfo.MaxAbsRowSum = NaN;

    case "power05"
        XmodelAll = local_powerTransformAbundance(Xraw, 0.5);
        abundanceInfo = struct();
        abundanceInfo.Transform = "power05";
        abundanceInfo.InputScale = "proportions";
        abundanceInfo.Pseudocount = 0;
        abundanceInfo.MaxAbsRowSum = NaN;
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

idxTest = false(size(Y));

for iType = 1:numel(communityTypes)
    idxType = find(clustID == communityTypes(iType));
    nTestType = max(1, round(holdoutFraction * numel(idxType)));
    idxType = idxType(randperm(numel(idxType)));
    idxTest(idxType(1:nTestType)) = true;
end

idxTrain = ~idxTest;
dataSet = strings(size(Y));
dataSet(idxTrain) = "train";
dataSet(idxTest) = "test";

sampleWeights = ones(size(Y));

if useClusterBalancedRidgeWeights
    sampleWeights(:) = 0;

    for iType = 1:numel(communityTypes)
        idxTypeTrain = clustID == communityTypes(iType) & idxTrain;

        if any(idxTypeTrain)
            sampleWeights(idxTypeTrain) = 1 / sum(idxTypeTrain);
        end
    end

    sampleWeights(idxTrain) = sampleWeights(idxTrain) ./ mean(sampleWeights(idxTrain), 'omitnan');
    sampleWeights(idxTest) = 1;
end

idxVariableSpecies = std(X(idxTrain, :), 0, 1, 'omitnan') > 0;
X = X(:, idxVariableSpecies);

modelPredictorNames = predictorNames(idxVariableSpecies);
modelSpeciesNames = speciesNames(idxVariableSpecies);
modelSpeciesPhyla = speciesPhyla(idxVariableSpecies);
modelSpeciesRoles = speciesRoles(idxVariableSpecies);
modelSpeciesPlotNames = speciesPlotNames(idxVariableSpecies);

xMu = mean(X(idxTrain, :), 1, 'omitnan');
xSigma = std(X(idxTrain, :), 0, 1, 'omitnan');
xSigma(xSigma == 0 | isnan(xSigma)) = 1;
Xz = (X - xMu) ./ xSigma;

nObs = numel(Y);
nPredictors = size(Xz, 2);
nModelFolds = min(nFolds, sum(idxTrain));
cvMdl = [];
cvMSE = [];
cv_inner = [];
bestLambda = NaN;
bestCVMSE = NaN;
cvRMSETransformed = NaN;
betaStandardized = zeros(nPredictors, 1);
biasStandardized = 0;
isSelectedPredictor = false(size(modelPredictorNames));

switch modelMethod
    case "stepwise"
        nModelFolds = NaN;

        tStepwiseTrain = array2table(Xz(idxTrain, :), 'VariableNames', cellstr(modelPredictorNames));
        tStepwiseTrain.ChlorophyllTransformed = Y(idxTrain);
        tStepwiseAll = array2table(Xz, 'VariableNames', cellstr(modelPredictorNames));

        linearMdl = stepwiselm(tStepwiseTrain, 'constant', ...
            'Upper', 'linear', ...
            'ResponseVar', 'ChlorophyllTransformed', ...
            'PredictorVars', cellstr(modelPredictorNames), ...
            'Criterion', char(stepwiseCriterion), ...
            'PEnter', stepwisePEnter, ...
            'PRemove', stepwisePRemove, ...
            'Verbose', 0);

        YhatTransformed = predict(linearMdl, tStepwiseAll);

        coefTbl = linearMdl.Coefficients;
        coefNames = string(coefTbl.Properties.RowNames);
        idxIntercept = coefNames == "(Intercept)" | lower(coefNames) == "intercept";

        if any(idxIntercept)
            biasStandardized = coefTbl.Estimate(find(idxIntercept, 1));
        end

        [isSelectedPredictor, idxCoef] = ismember(modelPredictorNames, coefNames);
        betaStandardized(isSelectedPredictor) = coefTbl.Estimate(idxCoef(isSelectedPredictor));

    case "ridge"
        if blockedCV
            nModelFolds = min(nModelFolds, numel(unique(cruiseID(idxTrain))));
            [~, ~, ~, ~, cv_inner] = ml_blockedCV(Xz(idxTrain, :), Y(idxTrain), ...
                cruiseID(idxTrain), monthNum(idxTrain), nModelFolds);
            nModelFolds = cv_inner.NumTestSets;

            cvMdl = fitrlinear(Xz(idxTrain, :), Y(idxTrain), ...
                'Learner', 'leastsquares', ...
                'Regularization', 'ridge', ...
                'CVPartition', cv_inner, ...
                'Lambda', lambdaGrid, ...
                'Weights', sampleWeights(idxTrain));
        else
            cvMdl = fitrlinear(Xz(idxTrain, :), Y(idxTrain), ...
                'Learner', 'leastsquares', ...
                'Regularization', 'ridge', ...
                'KFold', nModelFolds, ...
                'Lambda', lambdaGrid, ...
                'Weights', sampleWeights(idxTrain));
        end

        cvMSE = kfoldLoss(cvMdl);
        minCVMSE = min(cvMSE);
        idxAllowedLambda = cvMSE <= minCVMSE * (1 + lambdaSelectionTolerance);
        idxBestLambda = find(idxAllowedLambda, 1, 'last');
        bestCVMSE = cvMSE(idxBestLambda);
        bestLambda = lambdaGrid(idxBestLambda);
        cvRMSETransformed = sqrt(bestCVMSE);

        linearMdl = fitrlinear(Xz(idxTrain, :), Y(idxTrain), ...
            'Learner', 'leastsquares', ...
            'Regularization', 'ridge', ...
            'Lambda', bestLambda, ...
            'Weights', sampleWeights(idxTrain));

        YhatTransformed = predict(linearMdl, Xz);
        betaStandardized = linearMdl.Beta(:);
        biasStandardized = linearMdl.Bias;
        isSelectedPredictor = abs(betaStandardized) > 0;

    otherwise
        error('Unknown modelMethod: %s.', modelMethod);
end

YhatOriginal = chlBackTransform(YhatTransformed);

rSquaredTrainTransformed = hf_stat_R2(Y(idxTrain), YhatTransformed(idxTrain));
rSquaredTestTransformed = hf_stat_R2(Y(idxTest), YhatTransformed(idxTest));
rSquaredTrainOriginal = hf_stat_R2(Yoriginal(idxTrain), YhatOriginal(idxTrain));
rSquaredTestOriginal = hf_stat_R2(Yoriginal(idxTest), YhatOriginal(idxTest));
rmseTrainTransformed = sqrt(mean((Y(idxTrain) - YhatTransformed(idxTrain)).^2, 'omitnan'));
rmseTestTransformed = sqrt(mean((Y(idxTest) - YhatTransformed(idxTest)).^2, 'omitnan'));
rmseTrainOriginal = sqrt(mean((Yoriginal(idxTrain) - YhatOriginal(idxTrain)).^2, 'omitnan'));
rmseTestOriginal = sqrt(mean((Yoriginal(idxTest) - YhatOriginal(idxTest)).^2, 'omitnan'));

%% Contributions by community type

nSelectedPredictors = nnz(isSelectedPredictor);
betaUnstandardized = betaStandardized ./ xSigma(:);
biasUnstandardized = biasStandardized - sum(betaStandardized(:) .* xMu(:) ./ xSigma(:), 'omitnan');

XrawModel = Xraw(validRows, idxVariableSpecies);
Xunstandardized = X;
sampleContrib = Xunstandardized .* betaUnstandardized';

medianContribution = nan(numel(modelSpeciesNames), numel(communityTypes));
contributionQLow = nan(numel(modelSpeciesNames), numel(communityTypes));
contributionQHigh = nan(numel(modelSpeciesNames), numel(communityTypes));
meanRelAbund = nan(numel(modelSpeciesNames), numel(communityTypes));
meanUnstandardizedPredictor = nan(numel(modelSpeciesNames), numel(communityTypes));
meanStandardizedPredictor = nan(numel(modelSpeciesNames), numel(communityTypes));
nSamplesType = nan(numel(communityTypes), 1);

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    idxType = clustID == typeID;
    nSamplesType(iType) = sum(idxType);

    meanRelAbund(:, iType) = mean(XrawModel(idxType, :), 1, 'omitnan')';
    meanUnstandardizedPredictor(:, iType) = mean(Xunstandardized(idxType, :), 1, 'omitnan')';
    meanStandardizedPredictor(:, iType) = mean(Xz(idxType, :), 1, 'omitnan')';

    for iSpecies = 1:numel(modelSpeciesNames)
        vals = sampleContrib(idxType, iSpecies);
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

tCoef = table();
tCoef.Phylum = modelSpeciesPhyla;
tCoef.Species = modelSpeciesNames;
tCoef.Role = modelSpeciesRoles;
tCoef.PredictorName = modelPredictorNames;
tCoef.BetaStandardized = betaStandardized;
tCoef.BetaUnstandardized = betaUnstandardized;
tCoef.BiasStandardized = repmat(biasStandardized, nPredictors, 1);
tCoef.BiasUnstandardized = repmat(biasUnstandardized, nPredictors, 1);
tCoef.XMean = xMu(:);
tCoef.XSigma = xSigma(:);
tCoef.ModelMethod = repmat(modelMethod, nPredictors, 1);
tCoef.StepwiseSelected = isSelectedPredictor(:);
tCoef.StepwiseCriterion = repmat(stepwiseCriterion, nPredictors, 1);
tCoef.StepwisePEnter = repmat(stepwisePEnter, nPredictors, 1);
tCoef.StepwisePRemove = repmat(stepwisePRemove, nPredictors, 1);
tCoef.BestLambda = repmat(bestLambda, nPredictors, 1);

idxOutputSpecies = isSelectedPredictor(:);

if ~any(idxOutputSpecies)
    warning('The selected model did not retain any species predictors.');
end

tCoef = tCoef(idxOutputSpecies, :);
modelSpeciesNames = modelSpeciesNames(idxOutputSpecies);
modelSpeciesPhyla = modelSpeciesPhyla(idxOutputSpecies);
modelSpeciesRoles = modelSpeciesRoles(idxOutputSpecies);
modelSpeciesPlotNames = modelSpeciesPlotNames(idxOutputSpecies);
medianContribution = medianContribution(idxOutputSpecies, :);
contributionQLow = contributionQLow(idxOutputSpecies, :);
contributionQHigh = contributionQHigh(idxOutputSpecies, :);
meanRelAbund = meanRelAbund(idxOutputSpecies, :);
meanUnstandardizedPredictor = meanUnstandardizedPredictor(idxOutputSpecies, :);
meanStandardizedPredictor = meanStandardizedPredictor(idxOutputSpecies, :);

tContribution = tCoef(:, {'Phylum', 'Species', 'Role', 'PredictorName', ...
    'BetaStandardized', 'BetaUnstandardized', 'ModelMethod', ...
    'StepwiseSelected', 'StepwiseCriterion', 'StepwisePEnter', 'StepwisePRemove', ...
    'BestLambda'});

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    tContribution.(['MeanRelAbund_Type' num2str(typeID)]) = meanRelAbund(:, iType);
    tContribution.(['MeanUnstandardizedPredictor_Type' num2str(typeID)]) = meanUnstandardizedPredictor(:, iType);
    tContribution.(['MeanStandardizedPredictor_Type' num2str(typeID)]) = meanStandardizedPredictor(:, iType);
    tContribution.(['MedianContribution_Type' num2str(typeID)]) = medianContribution(:, iType);
    tContribution.(['Contribution' char(contributionIntervalLabel) 'Low_Type' num2str(typeID)]) = contributionQLow(:, iType);
    tContribution.(['Contribution' char(contributionIntervalLabel) 'High_Type' num2str(typeID)]) = contributionQHigh(:, iType);
    tContribution.(['NSamples_Type' num2str(typeID)]) = repmat(nSamplesType(iType), height(tContribution), 1);
end

tModelSummary = table(nObs, sum(idxTrain), sum(idxTest), holdoutFraction, ...
    nPredictors, nSelectedPredictors, nModelFolds, ...
    modelMethod, stepwiseCriterion, stepwisePEnter, stepwisePRemove, ...
    blockedCV, useClusterBalancedRidgeWeights, chlorophyllTransform, abundanceInfo.Transform, abundanceInfo.InputScale, abundanceInfo.Pseudocount, ...
    abundanceInfo.MaxAbsRowSum, lambdaSelectionTolerance, bestLambda, cvRMSETransformed, ...
    rSquaredTrainTransformed, rSquaredTestTransformed, rmseTrainTransformed, rmseTestTransformed, ...
    rSquaredTrainOriginal, rSquaredTestOriginal, rmseTrainOriginal, rmseTestOriginal, ...
    'VariableNames', {'NSamples', 'NTrain', 'NTest', 'HoldoutFraction', ...
    'NPredictors', 'NSelectedPredictors', 'KFolds', ...
    'ModelMethod', 'StepwiseCriterion', 'StepwisePEnter', 'StepwisePRemove', ...
    'BlockedCV', 'UseClusterBalancedRidgeWeights', 'ChlorophyllTransform', 'AbundanceTransform', ...
    'AbundanceInputScale', 'AbundancePseudoCount', 'CLRMaxAbsRowSum', ...
    'LambdaSelectionTolerance', 'BestLambda', 'CV_RMSE_Transformed', ...
    'TrainR2_Transformed', 'TestR2_Transformed', 'TrainRMSE_Transformed', 'TestRMSE_Transformed', ...
    'TrainR2_Original', 'TestR2_Original', 'TrainRMSE_Original', 'TestRMSE_Original'});

tPred = table();
tPred.SampleIndex = find(validRows);
tPred.ClustID = clustID;
tPred.ChlorophyllObserved = Yoriginal;
tPred.ChlorophyllPredicted = YhatOriginal;
tPred.ChlorophyllObservedTransformed = Y;
tPred.ChlorophyllPredictedTransformed = YhatTransformed;
tPred.DataSet = dataSet;
tPred.SampleWeight = sampleWeights;

%% Save results

modelMethodLabel = char(upper(extractBefore(modelMethod + "_", 2)) + extractAfter(modelMethod, 1));
outModelFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllGlobal' modelMethodLabel 'ContributionModel.mat']);
outSummaryFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllGlobal' modelMethodLabel 'ContributionSummary.csv']);
outCoefFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllGlobal' modelMethodLabel 'ContributionCoefficients.csv']);
outContributionFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllGlobal' modelMethodLabel 'ContributionByType.csv']);
outPredFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllGlobal' modelMethodLabel 'ContributionPredictions.csv']);

save(outModelFile, 'linearMdl', 'cvMdl', 'tModelSummary', 'tCoef', ...
    'tContribution', 'tPred', 'coreSpeciesData', 'chlorophyllTransform', ...
    'modelMethod', 'stepwiseCriterion', 'stepwisePEnter', 'stepwisePRemove', ...
    'lambdaGrid', 'lambdaSelectionTolerance', 'blockedCV', 'cv_inner', ...
    'useClusterBalancedRidgeWeights', 'sampleWeights', 'cvMSE', 'bestLambda', 'bestCVMSE', ...
    'contributionCILevel', 'contributionQuantileProbs', 'contributionIntervalLabel', ...
    'xMu', 'xSigma', 'idxOutputSpecies', ...
    'modelAbundanceTransform', 'abundanceInfo', 'medianContribution', ...
    'contributionQLow', 'contributionQHigh');
writetable(tModelSummary, outSummaryFile);
writetable(tCoef, outCoefFile);
writetable(tContribution, outContributionFile);
writetable(tPred, outPredFile);

%% Plot ridge cross-validation loss

if modelMethod == "ridge" && ~isempty(cvMSE)
    [fgLoss, MaxFigure] = addFindFigure('Global ridge lambda CV loss', ...
        MaxFigure, [40, 80, 720, 470], 'normal', 'on');
    clf(fgLoss);
    fgLoss.Color = 'w';

    axLoss = axes(fgLoss);
    semilogx(axLoss, lambdaGrid, cvMSE, '-o', ...
        'Color', [0.20 0.38 0.80], ...
        'MarkerFaceColor', [0.20 0.38 0.80], ...
        'MarkerSize', 4, ...
        'LineWidth', 1.2);
    hold(axLoss, 'on');
    xline(axLoss, bestLambda, 'k--', ...
        sprintf('\\lambda = %.3g', bestLambda), ...
        'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom');
    grid(axLoss, 'on');
    xlabel(axLoss, '\lambda');
    ylabel(axLoss, '5-fold CV MSE');
    title(axLoss, sprintf('Global ridge CV loss, best \\lambda = %.3g', bestLambda));
end

%% Plot contribution estimates by community type

[fgB, MaxFigure] = addFindFigure(['Global ' char(modelMethod) ' core species contributions'], ...
    MaxFigure, [24 332 1461 455], 'normal', 'on');

barOpts = struct();
barOpts.Figure = fgB;
barOpts.ClearFigure = true;
barOpts.XLabel = 'Community type';
barOpts.YLabel = char(coreSpeciesData.GroupName);
barOpts.Title = sprintf('Global %s species contribution to chlorophyll (linear scale, median and %.0f-%.0f%% quantiles)', ...
    char(modelMethod), contributionQuantileProbs(1) * 100, contributionQuantileProbs(3) * 100);
barOpts.ColorbarLabel = 'median unstandardized predictor \times beta (linear scale)';
barOpts.RoleColors = roleClrs;
barOpts.RoleShortNames = roleShortNames;
barOpts.RoleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};

[fgB, axB, cb, maxAbsPlotContribution] = plotGroupedIntervalBarsHor( ...
    medianContribution, contributionQLow, contributionQHigh, ...
    compose('Type %d', communityTypes), modelSpeciesPlotNames, ...
    modelSpeciesPhyla, modelSpeciesRoles, barOpts);

%% Observed/predicted diagnostics

for iSet = 1:2
    plotSet = ["train", "test"];
    plotSet = plotSet(iSet);
    idxSet = tPred.DataSet == plotSet;

    [fgPred, MaxFigure] = addFindFigure(['Global ' char(modelMethod) ' core species observed predicted ' char(plotSet)], ...
        MaxFigure, [10, 20, 1100 650], 'normal', 'on');
    clf(fgPred);
    fgPred.Color = 'w';

    tlPred = tiledlayout(fgPred, 'flow', 'Padding', 'compact', 'TileSpacing', 'compact');
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
        title(axPred, sprintf('%s Type %d, R^2_{log}=%.2f, r_{log}=%.2f', ...
            char(plotSet), typeID, r2Log, rLog));
    end
end

%% Mean contribution across community types

if ~isempty(modelSpeciesNames)
    contributionByType = medianContribution;
    meanContributionAcrossTypes = mean(contributionByType, 2, 'omitnan');
    nContributionTypes = sum(~isnan(contributionByType), 2);
    semContributionAcrossTypes = std(contributionByType, 0, 2, 'omitnan') ./ sqrt(nContributionTypes);
    semContributionAcrossTypes(nContributionTypes < 2 | isnan(semContributionAcrossTypes)) = 0;
    ciHalfWidth = 1.96 * semContributionAcrossTypes;
    contributionCILow = meanContributionAcrossTypes - ciHalfWidth;
    contributionCIHigh = meanContributionAcrossTypes + ciHalfWidth;

    [~, idxMeanContributionSort] = sort(meanContributionAcrossTypes, 'descend', 'MissingPlacement', 'last');

    sortedMeanContribution = meanContributionAcrossTypes(idxMeanContributionSort);
    sortedContributionLow = contributionCILow(idxMeanContributionSort);
    sortedContributionHigh = contributionCIHigh(idxMeanContributionSort);
    sortedSpeciesNames = modelSpeciesPlotNames(idxMeanContributionSort);
    sortedSpeciesRoles = modelSpeciesRoles(idxMeanContributionSort);

    maxAbsMeanContribution = max(abs([sortedContributionLow(:); sortedContributionHigh(:)]), [], 'omitnan');

    if isempty(maxAbsMeanContribution) || maxAbsMeanContribution == 0 || isnan(maxAbsMeanContribution)
        maxAbsMeanContribution = 1;
    end

    cmap = local_blueWhiteRed(256);

    [fgMean, MaxFigure] = addFindFigure(['Global ' char(modelMethod) ' mean species contribution'], ...
        MaxFigure, [60 120 900 720], 'normal', 'on');
    clf(fgMean);
    fgMean.Color = 'w';

    axMean = axes(fgMean);
    hold(axMean, 'on');
    axMean.Position = [0.32 0.08 0.56 0.84];

    yMean = (1:numel(sortedMeanContribution))';
    barColors = nan(numel(sortedMeanContribution), 3);

    for iSpecies = 1:numel(sortedMeanContribution)
        barColors(iSpecies, :) = local_valueColor(sortedMeanContribution(iSpecies), maxAbsMeanContribution, cmap);
    end

    bMean = barh(axMean, yMean, sortedMeanContribution, 0.62, ...
        'FaceColor', 'flat', ...
        'EdgeColor', 'none');
    bMean.CData = barColors;

    for iSpecies = 1:numel(sortedMeanContribution)
        plot(axMean, [sortedContributionLow(iSpecies) sortedContributionHigh(iSpecies)], ...
            [yMean(iSpecies) yMean(iSpecies)], ...
            'k-', 'LineWidth', 1.0);
        plot(axMean, [sortedContributionLow(iSpecies) sortedContributionLow(iSpecies)], ...
            yMean(iSpecies) + [-0.13 0.13], ...
            'k-', 'LineWidth', 1.0);
        plot(axMean, [sortedContributionHigh(iSpecies) sortedContributionHigh(iSpecies)], ...
            yMean(iSpecies) + [-0.13 0.13], ...
            'k-', 'LineWidth', 1.0);
    end

    xline(axMean, 0, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
    set(axMean, 'YDir', 'reverse');
    yticks(axMean, yMean);
    yticklabels(axMean, local_coloredYTickLabels(sortedSpeciesNames, sortedSpeciesRoles, roleClrs));
    axMean.TickLabelInterpreter = 'tex';
    axMean.FontSize = 8;
    ylim(axMean, [0.5 numel(sortedMeanContribution) + 0.5]);
    xlim(axMean, [-maxAbsMeanContribution maxAbsMeanContribution] * 1.08);
    xlabel(axMean, 'Mean contribution across community types');
    ylabel(axMean, char(coreSpeciesData.GroupName));
    title(axMean, sprintf('Global %s mean species contribution across types (mean \\pm 95%% CI)', char(modelMethod)));
    grid(axMean, 'on');
    box(axMean, 'on');
    colormap(axMean, cmap);
    clim(axMean, [-maxAbsMeanContribution maxAbsMeanContribution]);
    cbMean = colorbar(axMean);
    cbMean.Label.String = 'mean contribution';
    cbMean.Position = [0.90 0.17 0.018 0.66];

    local_addRoleLegend(axMean, roleClrs, roleShortNames, 8);
end

disp('Finished.');

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

    info = struct();
    info.Transform = "clr";
    info.InputScale = "proportions";
    info.Pseudocount = pseudocount;
    info.MaxAbsRowSum = max(abs(sum(Xclr, 2, 'omitnan')), [], 'omitnan');
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

    for iArg = 1:nargin
        x = varargin{iArg};
        vals = [vals; x(:)]; %#ok<AGROW>
    end

    vals = vals(vals > 0 & isfinite(vals));

    if isempty(vals)
        lim = [1e-3 1];
        return
    end

    mn = min(vals);
    mx = max(vals);

    if mx <= mn
        mn = mn / 1.2;
        mx = mx * 1.2;
    else
        pad = 0.06 * (log10(mx) - log10(mn));
        mn = 10 .^ (log10(mn) - pad);
        mx = 10 .^ (log10(mx) + pad);
    end

    lim = [mn mx];
end

function clr = local_roleColor(roleName, roleClrs)

    roleName = char(string(roleName));

    if isKey(roleClrs, roleName)
        clr = roleClrs(roleName);
    else
        clr = roleClrs('Unknown');
    end
end

function tickLabels = local_coloredYTickLabels(labels, roleNames, roleClrs)

    labels = string(labels);
    roleNames = string(roleNames);
    tickLabels = cell(numel(labels), 1);

    for i = 1:numel(labels)
        clr = local_roleColor(roleNames(i), roleClrs);
        label = local_escapeTex(labels(i));

        tickLabels{i} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%s', ...
            clr(1), clr(2), clr(3), label);
    end
end

function label = local_escapeTex(label)

    label = char(string(label));
    label = strrep(label, '\', '\\');
    label = strrep(label, '_', '\_');
    label = strrep(label, '^', '\^');
    label = strrep(label, '{', '\{');
    label = strrep(label, '}', '\}');
end

function local_addRoleLegend(ax, roleClrs, roleShortNames, fontSize)

    roleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};
    x0 = -0.13;
    y0 = 1.055;
    dx = 0.075;

    for i = 1:numel(roleOrder)
        roleName = roleOrder{i};
        clr = roleClrs(roleName);
        label = roleShortNames(roleName);
        x = x0 + (i - 1) * dx;

        text(ax, x, y0, label, ...
            'Units', 'normalized', ...
            'Color', clr, ...
            'FontSize', fontSize, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
    end
end

function cmap = local_blueWhiteRed(n)

    if nargin < 1
        n = 256;
    end

    n1 = floor(n / 2);
    n2 = n - n1;

    blue = [0.20 0.38 0.80];
    white = [1.00 1.00 1.00];
    red = [0.82 0.18 0.12];

    cmap1 = [linspace(blue(1), white(1), n1)', ...
             linspace(blue(2), white(2), n1)', ...
             linspace(blue(3), white(3), n1)'];

    cmap2 = [linspace(white(1), red(1), n2)', ...
             linspace(white(2), red(2), n2)', ...
             linspace(white(3), red(3), n2)'];

    cmap = [cmap1; cmap2];
end

function clr = local_valueColor(value, maxAbsValue, cmap)

    if maxAbsValue <= 0 || isnan(maxAbsValue)
        maxAbsValue = 1;
    end

    t = (value / maxAbsValue + 1) / 2;
    t = max(0, min(1, t));
    idx = 1 + round(t * (size(cmap, 1) - 1));
    clr = cmap(idx, :);
end
