% Fram_0277_CoreSpecContribution_modelSweep.m
%
% Benchmark per-community-type chlorophyll models across model method,
% abundance transform, and chlorophyll Box-Cox transform.
% Saves only the final performance summary table.

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

%% Parameters

speciesNames = string(coreSpeciesData.Species(:));
predictorNamesAll = matlab.lang.makeValidName(speciesNames);
predictorNamesAll = matlab.lang.makeUniqueStrings(predictorNamesAll);

Xraw = coreSpeciesData.Abundance;
YoriginalAll = tStPredictors.chlorophyll_a_0;
clustIDAll = tDMC.ClustID(:);
communityTypes = sort(unique(clustIDAll(:)))';

modelMethods = ["stepwise", "ridge"];
abundanceTransforms = ["power025", "power05", "none"];
boxcoxLambdaGrid = [0, 0.25, 0.5];

stepwiseCriterion = "bic";
stepwisePEnter = 0.05;
stepwisePRemove = 0.10;
lambdaGrid = logspace(-2, 2, 100);
nFolds = 5;
rng(1);

%% Sweep

tSweep = table();

for iMethod = 1:numel(modelMethods)
    modelMethod = modelMethods(iMethod);

    for iAbund = 1:numel(abundanceTransforms)
        modelAbundanceTransform = abundanceTransforms(iAbund);
        [XmodelAll, abundanceInfo] = local_transformAbundance(Xraw, modelAbundanceTransform);

        for iLambdaChl = 1:numel(boxcoxLambdaGrid)
            boxcoxLambdaChl = boxcoxLambdaGrid(iLambdaChl);
            chlTransform = @(y) local_boxcoxTransform(y, boxcoxLambdaChl);
            chlBackTransform = @(y) local_boxcoxInverse(y, boxcoxLambdaChl);

            validRows = ~isnan(YoriginalAll) & YoriginalAll > 0 & all(~isnan(XmodelAll), 2);
            X = XmodelAll(validRows, :);
            Yoriginal = YoriginalAll(validRows);
            Y = chlTransform(Yoriginal);
            clustID = clustIDAll(validRows);

            idxVariableSpecies = std(X, 0, 1, 'omitnan') > 0;
            X = X(:, idxVariableSpecies);
            predictorNames = predictorNamesAll(idxVariableSpecies);
            nPredictors = size(X, 2);

            nTypes = numel(communityTypes);
            r2TransformedByType = nan(1, nTypes);
            r2OriginalByType = nan(1, nTypes);
            rmseTransformedByType = nan(1, nTypes);
            rmseOriginalByType = nan(1, nTypes);
            nSamplesByType = nan(1, nTypes);
            nSelectedByType = nan(1, nTypes);
            bestLambdaByType = nan(1, nTypes);
            cvRMSETransformedByType = nan(1, nTypes);

            for iType = 1:nTypes
                typeID = communityTypes(iType);
                idxType = clustID == typeID;
                nSamplesByType(iType) = sum(idxType);

                Xtype = X(idxType, :);
                Ytype = Y(idxType);
                YoriginalType = Yoriginal(idxType);

                if numel(Ytype) < 3
                    continue
                end

                nOuterFolds = min(nFolds, numel(Ytype));
                cvOuter = cvpartition(numel(Ytype), 'KFold', nOuterFolds);
                YhatTransformed = nan(size(Ytype));
                foldNSelected = nan(nOuterFolds, 1);
                foldBestLambda = nan(nOuterFolds, 1);
                foldInnerRMSE = nan(nOuterFolds, 1);

                for iFold = 1:nOuterFolds
                    idxTrain = training(cvOuter, iFold);
                    idxTest = test(cvOuter, iFold);

                    Xtrain = Xtype(idxTrain, :);
                    Xtest = Xtype(idxTest, :);
                    Ytrain = Ytype(idxTrain);

                    idxVariableTrain = std(Xtrain, 0, 1, 'omitnan') > 0;

                    if nnz(idxVariableTrain) == 0 || numel(Ytrain) < 3
                        YhatTransformed(idxTest) = mean(Ytrain, 'omitnan');
                        foldNSelected(iFold) = 0;
                        continue
                    end

                    Xtrain = Xtrain(:, idxVariableTrain);
                    Xtest = Xtest(:, idxVariableTrain);
                    predictorNamesTrain = predictorNames(idxVariableTrain);

                    xMu = mean(Xtrain, 1, 'omitnan');
                    xSigma = std(Xtrain, 0, 1, 'omitnan');
                    xSigma(xSigma == 0 | isnan(xSigma)) = 1;
                    XtrainZ = (Xtrain - xMu) ./ xSigma;
                    XtestZ = (Xtest - xMu) ./ xSigma;

                    switch modelMethod
                        case "stepwise"
                            tTrain = array2table(XtrainZ, 'VariableNames', cellstr(predictorNamesTrain));
                            tTrain.ChlorophyllTransformed = Ytrain;

                            mdl = stepwiselm(tTrain, 'constant', ...
                                'Upper', 'linear', ...
                                'ResponseVar', 'ChlorophyllTransformed', ...
                                'PredictorVars', cellstr(predictorNamesTrain), ...
                                'Criterion', char(stepwiseCriterion), ...
                                'PEnter', stepwisePEnter, ...
                                'PRemove', stepwisePRemove, ...
                                'Verbose', 0);

                            tTest = array2table(XtestZ, 'VariableNames', cellstr(predictorNamesTrain));
                            YhatTransformed(idxTest) = predict(mdl, tTest);
                            coefNames = string(mdl.Coefficients.Properties.RowNames);
                            foldNSelected(iFold) = nnz(ismember(predictorNamesTrain, coefNames));

                        case "ridge"
                            nInnerFolds = min(nFolds, numel(Ytrain));

                            cvMdl = fitrlinear(XtrainZ, Ytrain, ...
                                'Learner', 'leastsquares', ...
                                'Regularization', 'ridge', ...
                                'KFold', nInnerFolds, ...
                                'Lambda', lambdaGrid);

                            cvMSE = kfoldLoss(cvMdl);
                            [bestCVMSE, idxBestLambda] = min(cvMSE);
                            bestLambda = lambdaGrid(idxBestLambda);
                            foldBestLambda(iFold) = bestLambda;
                            foldInnerRMSE(iFold) = sqrt(bestCVMSE);

                            mdl = fitrlinear(XtrainZ, Ytrain, ...
                                'Learner', 'leastsquares', ...
                                'Regularization', 'ridge', ...
                                'Lambda', bestLambda);

                            YhatTransformed(idxTest) = predict(mdl, XtestZ);
                            foldNSelected(iFold) = nnz(abs(mdl.Beta(:)) > 0);

                        otherwise
                            error('Unknown modelMethod: %s.', modelMethod);
                    end
                end

                YhatOriginal = chlBackTransform(YhatTransformed);

                r2TransformedByType(iType) = local_rSquared(Ytype, YhatTransformed);
                r2OriginalByType(iType) = local_rSquared(YoriginalType, YhatOriginal);
                rmseTransformedByType(iType) = sqrt(mean((Ytype - YhatTransformed).^2, 'omitnan'));
                rmseOriginalByType(iType) = sqrt(mean((YoriginalType - YhatOriginal).^2, 'omitnan'));
                nSelectedByType(iType) = mean(foldNSelected, 'omitnan');
                bestLambdaByType(iType) = median(foldBestLambda, 'omitnan');
                cvRMSETransformedByType(iType) = mean(foldInnerRMSE, 'omitnan');
            end

            tRow = table();
            tRow.ModelMethod = modelMethod;
            tRow.AbundanceTransform = modelAbundanceTransform;
            tRow.BoxCoxLambdaChl = boxcoxLambdaChl;
            tRow.NPredictors = nPredictors;
            tRow.StepwiseCriterion = stepwiseCriterion;
            tRow.StepwisePEnter = stepwisePEnter;
            tRow.StepwisePRemove = stepwisePRemove;
            tRow.AbundanceInputScale = abundanceInfo.InputScale;
            tRow.AbundancePseudoCount = abundanceInfo.Pseudocount;
            tRow.CLRMaxAbsRowSum = abundanceInfo.MaxAbsRowSum;
            tRow.EvaluationMode = "outer_cv_test";
            tRow.MeanTestR2_Transformed = mean(r2TransformedByType, 'omitnan');
            tRow.MinTestR2_Transformed = min(r2TransformedByType, [], 'omitnan');
            tRow.MedianTestR2_Transformed = median(r2TransformedByType, 'omitnan');
            tRow.MeanTestR2_Original = mean(r2OriginalByType, 'omitnan');
            tRow.MinTestR2_Original = min(r2OriginalByType, [], 'omitnan');
            tRow.MedianTestR2_Original = median(r2OriginalByType, 'omitnan');
            tRow.MeanTestRMSE_Transformed = mean(rmseTransformedByType, 'omitnan');
            tRow.MeanTestRMSE_Original = mean(rmseOriginalByType, 'omitnan');
            tRow.MeanNSelectedPredictors = mean(nSelectedByType, 'omitnan');
            tRow.MinNSelectedPredictors = min(nSelectedByType, [], 'omitnan');
            tRow.MeanBestLambda = mean(bestLambdaByType, 'omitnan');
            tRow.MedianBestLambda = median(bestLambdaByType, 'omitnan');
            tRow.MeanInnerCV_RMSE_Transformed = mean(cvRMSETransformedByType, 'omitnan');

            for iType = 1:nTypes
                typeID = communityTypes(iType);
                tRow.(['TestR2t_Type' num2str(typeID)]) = r2TransformedByType(iType);
                tRow.(['TestR2_Type' num2str(typeID)]) = r2OriginalByType(iType);
                tRow.(['N_Type' num2str(typeID)]) = nSamplesByType(iType);
                tRow.(['NSelected_Type' num2str(typeID)]) = nSelectedByType(iType);
                tRow.(['BestLambda_Type' num2str(typeID)]) = bestLambdaByType(iType);
            end

            tSweep = [tSweep; tRow]; %#ok<AGROW>

            fprintf('Finished %s / %s / lambdaChl=%g: mean test R2t=%.3f, min test R2t=%.3f\n', ...
                modelMethod, modelAbundanceTransform, boxcoxLambdaChl, ...
                tRow.MeanTestR2_Transformed, tRow.MinTestR2_Transformed);
        end
    end
end

%% Save and print ranking

tSweep = sortrows(tSweep, {'MeanTestR2_Transformed', 'MinTestR2_Transformed'}, {'descend', 'descend'});

outFile = fullfile(HGP.Data.Path2ResultsData, 'CoreSpeciesChlorophyllPerTypeModelSweep.csv');
writetable(tSweep, outFile);

disp('Top model sweep configurations:');
disp(tSweep(1:min(10, height(tSweep)), ...
    {'ModelMethod', 'AbundanceTransform', 'BoxCoxLambdaChl', ...
    'MeanTestR2_Transformed', 'MinTestR2_Transformed', 'MeanTestR2_Original', ...
    'MinTestR2_Original', 'MeanNSelectedPredictors'}));

fprintf('Saved model sweep table: %s\n', outFile);
disp('Finished.');

function [Xout, info] = local_transformAbundance(Xraw, transformName)

    switch string(transformName)
        case "none"
            Xout = Xraw;
            info = struct();
            info.Transform = "none";
            info.InputScale = "proportions";
            info.Pseudocount = 0;
            info.MaxAbsRowSum = NaN;

        case "power025"
            Xout = local_powerTransformAbundance(Xraw, 0.25);
            info = struct();
            info.Transform = "power025";
            info.InputScale = "proportions";
            info.Pseudocount = 0;
            info.MaxAbsRowSum = NaN;

        case "power05"
            Xout = local_powerTransformAbundance(Xraw, 0.5);
            info = struct();
            info.Transform = "power05";
            info.InputScale = "proportions";
            info.Pseudocount = 0;
            info.MaxAbsRowSum = NaN;

        case "clr"
            [Xout, info] = local_clrTransformAbundance(Xraw);

        otherwise
            error('Unknown abundance transform: %s.', transformName);
    end
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

    info = struct();
    info.Transform = "clr";
    info.InputScale = "proportions";
    info.Pseudocount = pseudocount;
    info.MaxAbsRowSum = max(abs(sum(Xclr, 2, 'omitnan')), [], 'omitnan');
end

function yTrans = local_boxcoxTransform(y, lambda)

    y = double(y);

    if abs(lambda) < eps
        yTrans = log10(y);
    else
        yTrans = (y .^ lambda - 1) ./ lambda;
    end
end

function y = local_boxcoxInverse(yTrans, lambda)

    yTrans = double(yTrans);

    if abs(lambda) < eps
        y = 10 .^ yTrans;
    else
        base = lambda .* yTrans + 1;
        base(base <= 0) = realmin;
        y = base .^ (1 ./ lambda);
    end
end

function r2 = local_rSquared(y, yHat)

    idx = ~isnan(y) & ~isnan(yHat);
    y = y(idx);
    yHat = yHat(idx);

    if numel(y) < 2
        r2 = NaN;
        return
    end

    ssTot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
    ssRes = sum((y - yHat).^2, 'omitnan');

    if ssTot <= 0
        r2 = NaN;
    else
        r2 = 1 - ssRes / ssTot;
    end
end
