% Fram_0277_CoreSpecContribution_detailed.m
%
% Fit separate chlorophyll models inside each community type.
% Then estimate the contribution of each core species as transformed abundance
% multiplied by the type-specific unstandardized regression beta.

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

idxVariableSpecies = std(X(idxTrain, :), 0, 1, 'omitnan') > 0;
X = X(:, idxVariableSpecies);

modelPredictorNames = predictorNames(idxVariableSpecies);
modelSpeciesNames = speciesNames(idxVariableSpecies);
modelSpeciesPhyla = speciesPhyla(idxVariableSpecies);
modelSpeciesRoles = speciesRoles(idxVariableSpecies);
modelSpeciesPlotNames = speciesPlotNames(idxVariableSpecies);

XrawModel = Xraw(validRows, idxVariableSpecies);
Xunstandardized = X;
Xz = nan(size(X));
nPredictors = size(X, 2);
nTypes = numel(communityTypes);

betaStandardizedMat = zeros(nPredictors, nTypes);
betaUnstandardizedMat = zeros(nPredictors, nTypes);
biasStandardizedVec = nan(nTypes, 1);
biasUnstandardizedVec = nan(nTypes, 1);
bestLambdaVec = nan(nTypes, 1);
bestCVMSEVec = nan(nTypes, 1);
cvRMSETransformedVec = nan(nTypes, 1);
nModelFoldsVec = nan(nTypes, 1);
nSamplesType = nan(nTypes, 1);
nSelectedPredictorsVec = nan(nTypes, 1);
isSelectedPredictorMat = false(nPredictors, nTypes);
xMuMat = nan(nPredictors, nTypes);
xSigmaMat = nan(nPredictors, nTypes);
YhatTransformed = nan(size(Y));
linearModels = cell(nTypes, 1);
cvModels = cell(nTypes, 1);
cvInnerByType = cell(nTypes, 1);
cvMSEByType = cell(nTypes, 1);

medianContribution = nan(nPredictors, nTypes);
contributionQLow = nan(nPredictors, nTypes);
contributionQHigh = nan(nPredictors, nTypes);
meanRelAbund = nan(nPredictors, nTypes);
meanUnstandardizedPredictor = nan(nPredictors, nTypes);
meanStandardizedPredictor = nan(nPredictors, nTypes);

tCoef = table();
tModelSummary = table();

for iType = 1:nTypes
    typeID = communityTypes(iType);
    idxType = clustID == typeID;
    nSamplesType(iType) = sum(idxType);
    idxTypeTrain = idxType & idxTrain;
    idxTypeTest = idxType & idxTest;

    Xtype = X(idxType, :);
    Ytype = Y(idxType);
    YoriginalType = Yoriginal(idxType);
    XtypeTrain = X(idxTypeTrain, :);
    YtypeTrain = Y(idxTypeTrain);
    idxVariableType = std(XtypeTrain, 0, 1, 'omitnan') > 0;

    xMuType = mean(XtypeTrain(:, idxVariableType), 1, 'omitnan');
    xSigmaType = std(XtypeTrain(:, idxVariableType), 0, 1, 'omitnan');
    xSigmaType(xSigmaType == 0 | isnan(xSigmaType)) = 1;
    XzType = (Xtype(:, idxVariableType) - xMuType) ./ xSigmaType;
    Xz(idxType, idxVariableType) = XzType;
    XzTypeTrain = (XtypeTrain(:, idxVariableType) - xMuType) ./ xSigmaType;

    predictorNamesType = modelPredictorNames(idxVariableType);
    betaStandardizedType = zeros(nnz(idxVariableType), 1);
    biasStandardizedType = 0;
    isSelectedType = false(nnz(idxVariableType), 1);
    bestLambda = NaN;
    bestCVMSE = NaN;
    cvRMSETransformed = NaN;
    nModelFolds = min(nFolds, nSamplesType(iType));
    cvMdl = [];
    cv_inner = [];
    cvMSE = [];

    switch modelMethod
        case "stepwise"
            nModelFolds = NaN;
            tTypeTrain = array2table(XzTypeTrain, 'VariableNames', cellstr(predictorNamesType));
            tTypeTrain.ChlorophyllTransformed = YtypeTrain;
            tTypeAll = array2table(XzType, 'VariableNames', cellstr(predictorNamesType));

            linearMdl = stepwiselm(tTypeTrain, 'constant', ...
                'Upper', 'linear', ...
                'ResponseVar', 'ChlorophyllTransformed', ...
                'PredictorVars', cellstr(predictorNamesType), ...
                'Criterion', char(stepwiseCriterion), ...
                'PEnter', stepwisePEnter, ...
                'PRemove', stepwisePRemove, ...
                'Verbose', 0);

            YhatType = predict(linearMdl, tTypeAll);
            coefTbl = linearMdl.Coefficients;
            coefNames = string(coefTbl.Properties.RowNames);
            idxIntercept = coefNames == "(Intercept)" | lower(coefNames) == "intercept";

            if any(idxIntercept)
                biasStandardizedType = coefTbl.Estimate(find(idxIntercept, 1));
            end

            [isSelectedType, idxCoef] = ismember(predictorNamesType, coefNames);
            betaStandardizedType(isSelectedType) = coefTbl.Estimate(idxCoef(isSelectedType));

        case "ridge"
            nModelFolds = min(nFolds, numel(YtypeTrain));

            nBlockedFolds = min(nModelFolds, numel(unique(cruiseID(idxTypeTrain))));

            if blockedCV && nBlockedFolds >= 2
                nModelFolds = nBlockedFolds;
                [~, ~, ~, ~, cv_inner] = ml_blockedCV(XzTypeTrain, YtypeTrain, ...
                    cruiseID(idxTypeTrain), monthNum(idxTypeTrain), nModelFolds);
                nModelFolds = cv_inner.NumTestSets;

                cvMdl = fitrlinear(XzTypeTrain, YtypeTrain, ...
                    'Learner', 'leastsquares', ...
                    'Regularization', 'ridge', ...
                    'CVPartition', cv_inner, ...
                    'Lambda', lambdaGrid);
            else
                cvMdl = fitrlinear(XzTypeTrain, YtypeTrain, ...
                    'Learner', 'leastsquares', ...
                    'Regularization', 'ridge', ...
                    'KFold', nModelFolds, ...
                    'Lambda', lambdaGrid);
            end

            cvMSE = kfoldLoss(cvMdl);
            minCVMSE = min(cvMSE);
            idxAllowedLambda = cvMSE <= minCVMSE * (1 + lambdaSelectionTolerance);
            idxBestLambda = find(idxAllowedLambda, 1, 'last');
            bestCVMSE = cvMSE(idxBestLambda);
            bestLambda = lambdaGrid(idxBestLambda);
            cvRMSETransformed = sqrt(bestCVMSE);

            linearMdl = fitrlinear(XzTypeTrain, YtypeTrain, ...
                'Learner', 'leastsquares', ...
                'Regularization', 'ridge', ...
                'Lambda', bestLambda);

            YhatType = predict(linearMdl, XzType);
            betaStandardizedType = linearMdl.Beta(:);
            biasStandardizedType = linearMdl.Bias;
            isSelectedType = abs(betaStandardizedType) > 0;

        otherwise
            error('Unknown modelMethod: %s.', modelMethod);
    end

    betaUnstandardizedType = betaStandardizedType ./ xSigmaType(:);
    biasUnstandardizedType = biasStandardizedType - ...
        sum(betaStandardizedType(:) .* xMuType(:) ./ xSigmaType(:), 'omitnan');

    betaStandardizedMat(idxVariableType, iType) = betaStandardizedType;
    betaUnstandardizedMat(idxVariableType, iType) = betaUnstandardizedType;
    isSelectedPredictorMat(idxVariableType, iType) = isSelectedType;
    xMuMat(idxVariableType, iType) = xMuType(:);
    xSigmaMat(idxVariableType, iType) = xSigmaType(:);
    biasStandardizedVec(iType) = biasStandardizedType;
    biasUnstandardizedVec(iType) = biasUnstandardizedType;
    bestLambdaVec(iType) = bestLambda;
    bestCVMSEVec(iType) = bestCVMSE;
    cvRMSETransformedVec(iType) = cvRMSETransformed;
    nModelFoldsVec(iType) = nModelFolds;
    nSelectedPredictorsVec(iType) = nnz(isSelectedType);
    YhatTransformed(idxType) = YhatType;
    linearModels{iType} = linearMdl;
    cvModels{iType} = cvMdl;
    cvInnerByType{iType} = cv_inner;
    cvMSEByType{iType} = cvMSE;

    sampleContribType = Xtype .* betaUnstandardizedMat(:, iType)';
    meanRelAbund(:, iType) = mean(XrawModel(idxType, :), 1, 'omitnan')';
    meanUnstandardizedPredictor(:, iType) = mean(Xtype, 1, 'omitnan')';
    meanStandardizedPredictor(idxVariableType, iType) = mean(XzType, 1, 'omitnan')';

    for iSpecies = 1:nPredictors
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

    coefType = table();
    coefType.ClustID = repmat(typeID, nPredictors, 1);
    coefType.Phylum = modelSpeciesPhyla;
    coefType.Species = modelSpeciesNames;
    coefType.Role = modelSpeciesRoles;
    coefType.PredictorName = modelPredictorNames;
    coefType.BetaStandardized = betaStandardizedMat(:, iType);
    coefType.BetaUnstandardized = betaUnstandardizedMat(:, iType);
    coefType.BiasStandardized = repmat(biasStandardizedType, nPredictors, 1);
    coefType.BiasUnstandardized = repmat(biasUnstandardizedType, nPredictors, 1);
    coefType.XMean = xMuMat(:, iType);
    coefType.XSigma = xSigmaMat(:, iType);
    coefType.ModelMethod = repmat(modelMethod, nPredictors, 1);
    coefType.StepwiseSelected = isSelectedPredictorMat(:, iType);
    coefType.StepwiseCriterion = repmat(stepwiseCriterion, nPredictors, 1);
    coefType.StepwisePEnter = repmat(stepwisePEnter, nPredictors, 1);
    coefType.StepwisePRemove = repmat(stepwisePRemove, nPredictors, 1);
    coefType.BestLambda = repmat(bestLambda, nPredictors, 1);
    tCoef = [tCoef; coefType(isSelectedPredictorMat(:, iType), :)]; %#ok<AGROW>

    YhatOriginalType = chlBackTransform(YhatType);
    idxTypeTrainLocal = idxTrain(idxType);
    idxTypeTestLocal = idxTest(idxType);

    summaryType = table(typeID, nSamplesType(iType), sum(idxTypeTrain), sum(idxTypeTest), holdoutFraction, ...
        nPredictors, nSelectedPredictorsVec(iType), ...
        nModelFolds, modelMethod, stepwiseCriterion, stepwisePEnter, stepwisePRemove, ...
        blockedCV, chlorophyllTransform, abundanceInfo.Transform, abundanceInfo.InputScale, abundanceInfo.Pseudocount, ...
        abundanceInfo.MaxAbsRowSum, lambdaSelectionTolerance, bestLambda, cvRMSETransformed, ...
        hf_stat_R2(Ytype(idxTypeTrainLocal), YhatType(idxTypeTrainLocal)), ...
        hf_stat_R2(Ytype(idxTypeTestLocal), YhatType(idxTypeTestLocal)), ...
        sqrt(mean((Ytype(idxTypeTrainLocal) - YhatType(idxTypeTrainLocal)).^2, 'omitnan')), ...
        sqrt(mean((Ytype(idxTypeTestLocal) - YhatType(idxTypeTestLocal)).^2, 'omitnan')), ...
        hf_stat_R2(YoriginalType(idxTypeTrainLocal), YhatOriginalType(idxTypeTrainLocal)), ...
        hf_stat_R2(YoriginalType(idxTypeTestLocal), YhatOriginalType(idxTypeTestLocal)), ...
        sqrt(mean((YoriginalType(idxTypeTrainLocal) - YhatOriginalType(idxTypeTrainLocal)).^2, 'omitnan')), ...
        sqrt(mean((YoriginalType(idxTypeTestLocal) - YhatOriginalType(idxTypeTestLocal)).^2, 'omitnan')), ...
        'VariableNames', {'ClustID', 'NSamples', 'NTrain', 'NTest', 'HoldoutFraction', ...
        'NPredictors', 'NSelectedPredictors', 'KFolds', ...
        'ModelMethod', 'StepwiseCriterion', 'StepwisePEnter', 'StepwisePRemove', ...
        'BlockedCV', 'ChlorophyllTransform', 'AbundanceTransform', ...
        'AbundanceInputScale', 'AbundancePseudoCount', 'CLRMaxAbsRowSum', ...
        'LambdaSelectionTolerance', 'BestLambda', 'CV_RMSE_Transformed', ...
        'TrainR2_Transformed', 'TestR2_Transformed', 'TrainRMSE_Transformed', 'TestRMSE_Transformed', ...
        'TrainR2_Original', 'TestR2_Original', 'TrainRMSE_Original', 'TestRMSE_Original'});
    tModelSummary = [tModelSummary; summaryType]; %#ok<AGROW>
end

YhatOriginal = chlBackTransform(YhatTransformed);

tPred = table();
tPred.SampleIndex = find(validRows);
tPred.ClustID = clustID;
tPred.ChlorophyllObserved = Yoriginal;
tPred.ChlorophyllPredicted = YhatOriginal;
tPred.ChlorophyllObservedTransformed = Y;
tPred.ChlorophyllPredictedTransformed = YhatTransformed;
tPred.DataSet = dataSet;

idxOutputSpecies = any(isSelectedPredictorMat, 2);

if ~any(idxOutputSpecies)
    warning('The selected models did not retain any species predictors.');
end

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

tContribution = table();
tContribution.Phylum = modelSpeciesPhyla;
tContribution.Species = modelSpeciesNames;
tContribution.Role = modelSpeciesRoles;
tContribution.PredictorName = modelPredictorNames(idxOutputSpecies);
tContribution.ModelMethod = repmat(modelMethod, nnz(idxOutputSpecies), 1);

for iType = 1:nTypes
    typeID = communityTypes(iType);
    tContribution.(['MeanRelAbund_Type' num2str(typeID)]) = meanRelAbund(:, iType);
    tContribution.(['MeanUnstandardizedPredictor_Type' num2str(typeID)]) = meanUnstandardizedPredictor(:, iType);
    tContribution.(['MeanStandardizedPredictor_Type' num2str(typeID)]) = meanStandardizedPredictor(:, iType);
    tContribution.(['MedianContribution_Type' num2str(typeID)]) = medianContribution(:, iType);
    tContribution.(['Contribution' char(contributionIntervalLabel) 'Low_Type' num2str(typeID)]) = contributionQLow(:, iType);
    tContribution.(['Contribution' char(contributionIntervalLabel) 'High_Type' num2str(typeID)]) = contributionQHigh(:, iType);
    tContribution.(['BetaUnstandardized_Type' num2str(typeID)]) = betaUnstandardizedMat(idxOutputSpecies, iType);
    tContribution.(['NSamples_Type' num2str(typeID)]) = repmat(nSamplesType(iType), height(tContribution), 1);
end

%% Save results

modelMethodLabel = char(upper(extractBefore(modelMethod + "_", 2)) + extractAfter(modelMethod, 1));
outModelFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllPerType' modelMethodLabel 'ContributionModel.mat']);
outSummaryFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllPerType' modelMethodLabel 'ContributionSummary.csv']);
outCoefFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllPerType' modelMethodLabel 'ContributionCoefficients.csv']);
outContributionFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllPerType' modelMethodLabel 'ContributionByType.csv']);
outPredFile = fullfile(HGP.Data.Path2ResultsData, ['CoreSpeciesChlorophyllPerType' modelMethodLabel 'ContributionPredictions.csv']);

save(outModelFile, 'linearModels', 'cvModels', 'cvInnerByType', 'cvMSEByType', 'tModelSummary', 'tCoef', ...
    'tContribution', 'tPred', 'coreSpeciesData', 'chlorophyllTransform', ...
    'modelMethod', 'stepwiseCriterion', 'stepwisePEnter', 'stepwisePRemove', ...
    'lambdaGrid', 'lambdaSelectionTolerance', 'blockedCV', 'bestLambdaVec', 'bestCVMSEVec', ...
    'contributionCILevel', 'contributionQuantileProbs', 'contributionIntervalLabel', ...
    'xMuMat', 'xSigmaMat', 'idxOutputSpecies', ...
    'betaStandardizedMat', 'betaUnstandardizedMat', 'biasStandardizedVec', 'biasUnstandardizedVec', ...
    'modelAbundanceTransform', 'abundanceInfo', 'medianContribution', ...
    'contributionQLow', 'contributionQHigh');
writetable(tModelSummary, outSummaryFile);
writetable(tCoef, outCoefFile);
writetable(tContribution, outContributionFile);
writetable(tPred, outPredFile);

%% Plot ridge cross-validation loss

if modelMethod == "ridge"
    [fgLoss, MaxFigure] = addFindFigure('Per-type ridge lambda CV loss', ...
        MaxFigure, [40, 80, 1100, 650], 'normal', 'on');
    clf(fgLoss);
    fgLoss.Color = 'w';

    tlLoss = tiledlayout(fgLoss, 'flow', ...
        'Padding', 'compact', ...
        'TileSpacing', 'compact');

    for iType = 1:numel(communityTypes)
        axLoss = nexttile(tlLoss);
        cvMSE = cvMSEByType{iType};

        if isempty(cvMSE)
            axis(axLoss, 'off');
            title(axLoss, sprintf('Type %d', communityTypes(iType)));
            continue
        end

        semilogx(axLoss, lambdaGrid, cvMSE, '-o', ...
            'Color', [0.20 0.38 0.80], ...
            'MarkerFaceColor', [0.20 0.38 0.80], ...
            'MarkerSize', 3, ...
            'LineWidth', 1.0);
        hold(axLoss, 'on');
        xline(axLoss, bestLambdaVec(iType), 'k--', ...
            sprintf('\\lambda = %.3g', bestLambdaVec(iType)), ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom');
        grid(axLoss, 'on');
        xlabel(axLoss, '\lambda');
        ylabel(axLoss, '5-fold CV MSE');
        title(axLoss, sprintf('Type %d ridge CV loss', communityTypes(iType)));
    end
end

%% Plot contribution estimates by community type

[fgB, MaxFigure] = addFindFigure(['Per-type ' char(modelMethod) ' core species contributions'], ...
    MaxFigure, [24 332 1461 455], 'normal', 'on');

barOpts = struct();
barOpts.Figure = fgB;
barOpts.ClearFigure = true;
barOpts.XLabel = 'Community type';
barOpts.YLabel = char(coreSpeciesData.GroupName);
barOpts.Title = sprintf('Per-type %s species contribution to chlorophyll (linear scale, median and %.0f-%.0f%% quantiles)', ...
    char(modelMethod), contributionQuantileProbs(1) * 100, contributionQuantileProbs(3) * 100);
barOpts.ColorbarLabel = 'median unstandardized predictor \times beta (linear scale)';
barOpts.RoleColors = roleClrs;
barOpts.RoleShortNames = roleShortNames;
barOpts.RoleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};

typeLabels = strings(size(communityTypes));
idxSetForLabels = tPred.DataSet == "test";

for iType = 1:numel(communityTypes)
    typeID = communityTypes(iType);
    idxType = tPred.ClustID == typeID & idxSetForLabels;
    idxValid = idxType & tPred.ChlorophyllObserved > 0 & tPred.ChlorophyllPredicted > 0;
    r2LogType = hf_stat_R2(log10(tPred.ChlorophyllObserved(idxValid)), log10(tPred.ChlorophyllPredicted(idxValid)));
    rCorr = corr(log10(tPred.ChlorophyllObserved(idxValid)), log10(tPred.ChlorophyllPredicted(idxValid)), 'Rows', 'complete');
    typeLabels(iType) = sprintf('Type %d \\newlineR^2_{log}=%.2f \\newline r_{log}=%.2f', typeID, r2LogType, rCorr);
end

[fgB, axB, cb, maxAbsPlotContribution] = plotGroupedIntervalBarsHor( ...
    medianContribution, contributionQLow, contributionQHigh, ...
    typeLabels, modelSpeciesPlotNames, ...
    modelSpeciesPhyla, modelSpeciesRoles, barOpts);

%% Observed/predicted diagnostics

plotSets = ["train", "test"];

for iSet = 1:numel(plotSets)
    plotSet = plotSets(iSet);
    idxSet = tPred.DataSet == plotSet;

    [fgPred, MaxFigure] = addFindFigure(['Per-type ' char(modelMethod) ' core species observed predicted ' char(plotSet)], ...
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

    [fgMean, MaxFigure] = addFindFigure(['Per-type ' char(modelMethod) ' mean species contribution'], ...
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
    title(axMean, sprintf('Per-type %s mean species contribution across types (mean \\pm 95%% CI)', char(modelMethod)));
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
