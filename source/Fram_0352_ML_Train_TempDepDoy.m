%Train the model for Temp, Doy and Depths to predict Cluster ID 

Fram_0350_ML_Init_TempDepDoy;

%Class model ID
ClassifID = 'ECOC_SVM';

trainParams = {};
classes = categories(Y_train);
K = numel(classes);

%%

switch ClassifID
    case 'ECOC_Auto'
        % X, Y: training data
        % cv_inner: cvpartition object (e.g. cvpartition(Y,'KFold',5))
        ml_train_func = @ml_train_ECOC;
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

        % Single output only
        MdlCV = fitcecoc( ...
            X_train, Y_train, ...
            'Learners',                tSVM, ...
            'Coding',                  'onevsone', ...
            'OptimizeHyperparameters', optimVars, ...
            'HyperparameterOptimizationOptions', optimOptions);

        % Get optimization results from the model property
        HyperparameterOptimizationResults = MdlCV.HyperparameterOptimizationResults;


        bestParams.boxconstr =HyperparameterOptimizationResults.bestPoint.BoxConstraint(1);
        bestParams.sigmas = HyperparameterOptimizationResults.bestPoint.KernelScale(1);
    case 'ECOC_Auto_Cost'

        Cost = ones(K) - eye(K);  % base 0/1 cost

        % Make errors for class 1 more expensive
        lambda = 1;  % try 2, 3, 5, etc.
        Cost(1, :) = lambda * Cost(1, :);

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

        % Single output only
        MdlCV = fitcecoc( ...
            X_train, Y_train, ...
            'Learners',                tSVM, ...
            'Coding',                  'onevsone', ...
            Cost_Options{:}, ...
            'OptimizeHyperparameters', optimVars, ...
            'HyperparameterOptimizationOptions', optimOptions);

        % Get optimization results from the model property
        HyperparameterOptimizationResults = MdlCV.HyperparameterOptimizationResults;


        bestParams.boxconstr =HyperparameterOptimizationResults.bestPoint.BoxConstraint(1);
        bestParams.sigmas = HyperparameterOptimizationResults.bestPoint.KernelScale(1);
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = logspace(log10(0.1), log10(30), 10);
        HP_grid.boxconstr  = logspace(-1, 2, 11);

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
    case 'RF'
        %% hyperparameter grid
        clearvars HP_grid
        HP_grid.NumTrees = [50, 100, 200, 300];
        HP_grid.MinLeafSize  = [5:5:20];
        P = numel(predictors);
        HP_grid.NumPredictorsToSample = unique([max(1, round(sqrt(P))), ceil(P/3), P]); % extern


        %handle of the training function ml_train_func(X, Y, params, Name, Value....)
        ml_train_func = @ml_train_RF;

        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train, Y_train, cv_inner, ml_train_func, HP_grid, trainParams{:});
        bestParams
end


%%

Mdl1 = ml_train_func(X_train(:, :), Y_train, bestParams,  'FitPosterior', true);

[cvErr, CVMdl] = mlTrainAndCVLoss(X_train(:, :), Y_train, cv_inner, ml_train_func, bestParams);
disp(['CV err ' NS(cvErr)])

% Suppose you trained cross-validated model, e.g.:

[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');

% 1) Out-of-fold predictions
yHat = kfoldPredict(CVMdl);   % each point predicted by a model that did not see it

% 2) Confusion matrix (counts)
C = confusionmat(Y_train, yHat);

% 3) Row-normalized confusion matrix (per-class mixing)
Cnorm = C ./ sum(C, 2);

% 4) Nice plot
confusionchart(Y_train, yHat, 'Normalization','row-normalized');  % mixing per true class

%Save the algorithm
save(fullfile(HGP.Data.Path2ResultsData, ['trainedModel_ClustIDonTempDoyDepth' LM_Model '.mat']), 'Mdl1');