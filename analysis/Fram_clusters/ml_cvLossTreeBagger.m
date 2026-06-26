function loss = ml_cvLossTreeBagger(params, X, Y, cv_blocked)
% params: table from bayesopt with fields:
%   NumTrees, MinLeafSize, NumPredictorsToSample
% X: predictors (matrix or table)
% Y: response vector
% cv_blocked: cvpartition object

    K = cv_blocked.NumTestSets;
    mseFold = zeros(K,1);

    for k = 1:K
        trainIdx = training(cv_blocked, k);
        testIdx  = test(cv_blocked, k);

        Xtr = X(trainIdx, :);
        Ytr = Y(trainIdx, :);
        Xte = X(testIdx, :);
        Yte = Y(testIdx, :);

        % Train random forest on this fold
        Mdl_fold = TreeBagger( ...
            params.NumTrees, Xtr, Ytr, ...
            'Method','regression', ...
            'MinLeafSize', params.MinLeafSize, ...
            'NumPredictorsToSample', params.NumPredictorsToSample, ...
            'OOBPrediction', 'off');   % on OOB in CV

        % Predict on validation fold
        Ypred = predict(Mdl_fold, Xte);    % numeric vector

        % Fold MSE
        mseFold(k) = mean((Yte - Ypred).^2, 'omitnan');
    end

    % Blocked-CV MSE для bayesopt
    loss = mean(mseFold, 'omitnan');
end
