function [Mdl, results] = ml_tune_model(LM_model,X_train,Y_train,cv_blocked)
switch LM_model
    case {'Boost1',  'Boost2'}
        TreeDepth  = find(strcmp(LM_model, {'Boost1',  'Boost2'}));
        maxNumSplits = 2^TreeDepth - 1;
        t = templateTree('MaxNumSplits', maxNumSplits);

        % Cross-validated ensemble using blocked CV
        [Mdl] = fitrensemble( ...
            X_train, Y_train, ...
            'Method','LSBoost', ...
            'Learners', t, ...
            'OptimizeHyperparameters', {'NumLearningCycles','LearnRate'}, ...
            'HyperparameterOptimizationOptions', struct( ...
            'CVPartition', cv_blocked, ...               % blocked CV
            'AcquisitionFunctionName','expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 50))  ;         % adjust as needed

        % Objective: k-fold loss (MSE by default)
        results.HPOptimization    = Mdl.HyperparameterOptimizationResults;
        results.bestCVLoss = results.HPOptimization.MinObjective;       % this is your CV loss
        results.bestPoint  = results.HPOptimization.bestPoint;    % best hyperparameters
    case 'RF'
            t = templateTree('Reproducible', true);
    
            [Mdl] = fitrensemble( ...
                X_train, Y_train, ...
                'Method', 'Bag', ...                 % bagging = random forest style
                'Learners', t, ...                   % tree learners
                'OptimizeHyperparameters', { ...
                'NumLearningCycles', ...         % number of trees
                'MinLeafSize', ...               % leaf size (controls depth/variance)
                'MaxNumSplits', ...              % max tree complexity
                'NumVariablesToSample' }, ...    % mtry
                'HyperparameterOptimizationOptions', struct( ...
                'CVPartition', cv_blocked, ...   % your blocked CV
                'AcquisitionFunctionName', 'expected-improvement-plus', ...
                'MaxObjectiveEvaluations', 50)); % adjust as needed

        % Hyperparameter optimization results (analog to your LSBoost code)
        results.HPOptimization = Mdl.HyperparameterOptimizationResults;
        results.bestCVLoss     = results.HPOptimization.MinObjective;       % CV loss
        results.bestPoint      = results.HPOptimization.XAtMinObjective;    % best hyperparams     

    case 'RFNative'
        numPred = size(X_train, 2);   % число предикторов (если X_train - таблица, будет работать так же)
        varsRF = [
            optimizableVariable('NumTrees',[50 500],'Type','integer')          %  NumLearningCycles
            optimizableVariable('MinLeafSize',[1 50],'Type','integer')         % Tree depth
            optimizableVariable('NumPredictorsToSample',[1 numPred],'Type','integer')  % mtry
            ];

        % --- 2. Objective-функция: blocked-CV MSE для TreeBagger ---

        objFun = @(params) ml_cvLossTreeBagger(params, X_train, Y_train, cv_blocked);

        % --- 3.optimization 

        resultsBO = bayesopt( ...
            objFun, varsRF, ...
            'IsObjectiveDeterministic', true, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 50, ...       % 
            'Verbose', 1, ...
            'PlotFcn', {});                          % 

        % --- 4. Save results 

        results.HPOptimization = resultsBO;
        results.bestCVLoss     = resultsBO.MinObjective;        %
        results.bestPoint      = resultsBO.XAtMinObjective;     %

        % --- 5.Final model 

        bp = results.bestPoint;

        Mdl = TreeBagger( ...
            bp.NumTrees, X_train, Y_train, ...
            'Method', 'regression', ...
            'MinLeafSize', bp.MinLeafSize, ...
            'NumPredictorsToSample', bp.NumPredictorsToSample, ...
            'OOBPrediction', 'on');     % just in case 
        
    case 'SVM_auto'
        Mdl = fitrsvm( X_train, Y_train, 'KernelFunction','gaussian','KernelScale','auto',...
            'Standardize',true);
        MdlCV = fitrsvm( X_train, Y_train, 'KernelFunction','gaussian', 'KernelScale',Mdl.KernelParameters.Scale, ...
            'Standardize',true, 'CVPartition', cv_blocked);
        results.bestCVLoss = kfoldLoss(MdlCV);
    case 'Lin_stepwise'
        Mdl = stepwiselm(X_train, Y_train, 'constant', ...
            'Upper', 'quadratic', ...     % quadratic + interactions
            'Criterion', 'BIC');                % silent
        modelFormula = char(Mdl.Formula);
        % 2) Cross-validated version of this stepwise model with blocked CV
        K = cv_blocked.NumTestSets;
        mseFold = zeros(K,1);

        for k = 1:K
            trainIdx = training(cv_blocked, k);
            testIdx  = test(cv_blocked, k);

            tblTrain = X_train(trainIdx, :);
            tblTrain.y = Y_train(trainIdx);
            tblTest  = X_train(testIdx, :);
            tblTest.y = Y_train(testIdx);

            % Fit linear model with fixed formula (no stepwise here)
            Mdl_fold = fitlm(tblTrain, modelFormula);

            % Predict on validation fold
            yPred = predict(Mdl_fold, tblTest);

            % Fold MSE
            mseFold(k) = mean((tblTest.y - yPred).^2, 'omitnan');
        end


        % 3) CV loss (MSE by default)
        results.bestCVLoss = mean(mseFold, 'omitnan');
    case 'FitNet'

        varsNN = [
            optimizableVariable('NumHidden',[5 20],'Type','integer')
            optimizableVariable('Lambda',[1e-6 1],'Transform','log')  % L2-regularization
            ];

        objFunNN = @(params) ml_cvLossNN(params, X_train, Y_train, cv_blocked);

        resultsBO = bayesopt( ...
            objFunNN, varsNN, ...
            'IsObjectiveDeterministic', true, ...
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 50, ...          % can change it
            'Verbose', 1, ...
            'PlotFcn', {});                             % without plots

        results.HPOptimization = resultsBO;
        results.bestCVLoss     = resultsBO.MinObjective;
        results.bestPoint      = resultsBO.XAtMinObjective;

        % final model
        bp = results.bestPoint;

       Mdl = fitrnet( ...
            X_train, Y_train, ...
            'LayerSizes', bp.NumHidden, ...
            'Lambda',     bp.Lambda, ...
            'Standardize', true);


    case 'GAM'
        % 1) Train GAM with hyperparameter optimization and blocked CV
        Mdl = fitrgam( ...
            X_train, Y_train, ...
            'OptimizeHyperparameters', 'auto', ...   % tune typical GAM hyperparams
            'HyperparameterOptimizationOptions', struct( ...
            'CVPartition', cv_blocked, ...       % blocked CV for optimization
            'AcquisitionFunctionName', 'expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 50, ...   % adjust as needed
            'ShowPlots', false, ...
            'Verbose', 1));

        % 2) Hyperparameter optimization results (аналогично fitrensemble)
        results.HPOptimization = Mdl.HyperparameterOptimizationResults;
        results.bestCVLoss     = results.HPOptimization.MinObjective;       % CV loss used in HPO
        results.bestPoint      = results.HPOptimization.XAtMinObjective;    % best hyperparameters


    case {'SVM'}
        [Mdl] = fitrsvm( ...
            X_train, Y_train, ...
            'Standardize', true, ...                  % standardize predictors
            'KernelFunction', 'gaussian', ...         % or 'linear', 'rbf' etc.
            'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ...
            'HyperparameterOptimizationOptions', struct( ...
            'CVPartition', cv_blocked, ...        % blocked CV
            'AcquisitionFunctionName','expected-improvement-plus', ...
            'MaxObjectiveEvaluations', 100));      % adjust as needed

        % Hyperparameter optimization results
        results.HPOptimization = Mdl.HyperparameterOptimizationResults;
        results.bestCVLoss     = results.HPOptimization.MinObjective;        % CV loss (objective)
        results.bestPoint      = results.HPOptimization.XAtMinObjective;

    otherwise
        error('LM_model not defined')
end