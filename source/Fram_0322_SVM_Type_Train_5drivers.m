%Feature selection

Fram_0310_SVM_Type_Init_5drivers;

%Class model ID
ClassifID = 'ECOC_SVM';

trainParams = {};
classes = categories(Y_train);
K = numel(classes);

%%

switch ClassifID
    case 'ECOC_SVM_Auto_Cost'
        GroupCounts = grpstats(table(Y_train, 'VariableNames', {'Y'}), 'Y', 'numel');
        GroupCounts.Weight = 1./GroupCounts.GroupCount;

        Cost = ones(K) - eye(K);
        for k = 1:K
            Cost(k, :) = Cost(k, :) * GroupCounts.Weight(k);
        end

        Cost_Options = {'Cost', Cost};

        %handle of the training function ml_train_func(X, Y, params, Name, Value....)
        ml_train_func = @(X, Y, params, varargin) ml_train_ECOC(X, Y, params, Cost_Options{:}, varargin{:});

        tSVM = templateSVM( ...
            'KernelFunction', 'rbf', ...
            'Standardize',    true, ...
            'Type',           'classification');

        optimVars = {'BoxConstraint', 'KernelScale'};

        optimOptions = struct( ...
            'Optimizer',               'bayesopt', ...
            'UseParallel',             true, ...
            'ShowPlots',               false, ...
            'Verbose',                 0, ...
            'CVPartition',             cv_inner, ...
            'MaxObjectiveEvaluations', 40);

        MdlCV = fitcecoc( ...
            X_train, Y_train, ...
            'Learners',                tSVM, ...
            'Coding',                  'onevsone', ...
            Cost_Options{:}, ...
            'OptimizeHyperparameters', optimVars, ...
            'HyperparameterOptimizationOptions', optimOptions);

        HyperparameterOptimizationResults = MdlCV.HyperparameterOptimizationResults;

        bestParams.boxconstr = HyperparameterOptimizationResults.bestPoint.BoxConstraint(1);
        bestParams.sigmas = HyperparameterOptimizationResults.bestPoint.KernelScale(1);
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = logspace(log10(1), log10(30), 20);
        HP_grid.boxconstr  = logspace(1, 2, 20);

        ix = @(c) find(strcmp(classes, string(c)));
        %
        K = numel(classes);
        GroupCounts = grpstats(table(Y_train, 'VariableNames', {'Y'}), 'Y', [ 'numel']);
        GroupCounts.Weight = 1./GroupCounts.GroupCount;


        Cost = ones(K) - eye(K);          % base 0/1 cost
        % define the cost
        for k = 1:K
            Cost(k, :) = Cost(k, :) * GroupCounts.Weight(k);
        end

        Cost_Options = {'Cost', Cost};

        %handle of the training function ml_train_func(X, Y, params, Name, Value....)
        ml_train_func = @ml_train_ECOC;

        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train, Y_train, cv_inner, ml_train_func, HP_grid);

        Sigmas   = paramMesh{1};   % MaxNumSplits
        BC  = paramMesh{2};   % NumLearningCycles
        [fg, MaxFigure] = addFindFigure('CV error, fit', MaxFigure, [10, 20, 400 400]);
        nexttile(1)
        contourf(BC, Sigmas, CVloss, 'ShowText', 'on'); colorbar;
        xlabel('BC')
        ylabel('Sigma')
        f_Set_scales({'logx', 'logy'})
        hold on
        scatter(bestParams.boxconstr, bestParams.sigmas, 30, 'red', 'filled');
        hold off
        drawnow

    otherwise
        error("ClassifID unknown")
end


%%

Mdl1 = ml_train_func(X_train(:, :), Y_train, bestParams,  'FitPosterior', true);

[cvErr, CVMdl] = mlTrainAndCVLoss(X_train(:, :), Y_train, cv_inner, ml_train_func, bestParams);
disp(['CV err ' NS(cvErr)])



% 1) Out-of-fold predictions from validation folds
yHat = kfoldPredict(CVMdl);   % each point predicted by a model that did not see it

% 2) Confusion matrix (counts)
C = confusionmat(Y_train, yHat, 'Order', classes);

% 3) Row-normalized confusion matrix (per-class mixing)
Cnorm = C ./ sum(C, 2);

% 4) Cross-validation summary statistics
classSupport = sum(C, 2);
classPredicted = sum(C, 1)';
classTP = diag(C);
classRecall = classTP ./ classSupport;
classPrecision = classTP ./ classPredicted;
classF1 = 2 .* classPrecision .* classRecall ./ (classPrecision + classRecall);

classRecall(~isfinite(classRecall)) = NaN;
classPrecision(~isfinite(classPrecision)) = NaN;
classF1(~isfinite(classF1)) = NaN;

CVStats = struct();
CVStats.ClassifID = ClassifID;
CVStats.BestParams = bestParams;
CVStats.CVPartition = cv_inner;
CVStats.CVLoss = cvErr;
foldLoss = kfoldLoss(CVMdl, 'Mode', 'individual');
CVStats.FoldLoss = foldLoss;
CVStats.CVLossMean = mean(foldLoss);
CVStats.CVLossStd = std(foldLoss);
CVStats.CVLossSE = std(foldLoss) ./ sqrt(numel(foldLoss));
CVStats.Accuracy = mean(yHat == Y_train);
CVStats.BalancedAccuracy = mean(classRecall, 'omitnan');
CVStats.ConfusionMatrix = C;
CVStats.ConfusionMatrixRowNormalized = Cnorm;
CVStats.ClassStats = table(classes(:), classSupport, classPredicted, ...
    classPrecision, classRecall, classF1, ...
    'VariableNames', {'Class', 'Support', 'Predicted', 'Precision', 'Recall', 'F1'});

disp(CVStats.ClassStats)
disp(['CV accuracy ' NS(CVStats.Accuracy)])
disp(['Cross-validation error was ' NS(CVStats.CVLossMean) ' +/- ' NS(CVStats.CVLossStd) ' SD across folds.'])
disp(['CV balanced accuracy ' NS(CVStats.BalancedAccuracy)])

Importance = permutationImportance(Mdl1);
predictorNamesUsed = renamePredictors(Importance.Predictor);

CVStats.PermutationImportance = Importance;
disp(Importance)

[~, impOrder] = sort(Importance.ImportanceMean, 'ascend');
impMean = Importance.ImportanceMean(impOrder);
impStd = Importance.ImportanceStandardDeviation(impOrder);

%% 5) Nice plot
[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1, 869 369], 'normal', 'on');
nexttile
cmChart = confusionchart(Y_train, yHat, 'Normalization','row-normalized');  % mixing per true class
predictorTitle = strjoin(predictorNamesUsed, ', ');
if strlength(predictorTitle) > 60
    splitIdx = ceil(numel(predictorNamesUsed)/2);
    predictorTitle = sprintf('%s\n%s', ...
        strjoin(predictorNamesUsed(1:splitIdx), ', '), ...
        strjoin(predictorNamesUsed(splitIdx+1:end), ', '));
end
cmChart.Title = sprintf('Predictors: %s', predictorTitle);
nexttile
%plot permutation feature importance


barh(impMean)
hold on
errorbar(impMean, 1:numel(impMean), impStd, 'horizontal', 'LineStyle', 'none', 'Color', 'k')
hold off
yticks(1:numel(impMean))
yticklabels(predictorNamesUsed(impOrder))
xlabel('Permutation importance')

AddLetters2Plots(fg, plotLabels, "HShift", -0.05)

%Save the algorithm
save(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']), 'Mdl1', 'CVStats');
