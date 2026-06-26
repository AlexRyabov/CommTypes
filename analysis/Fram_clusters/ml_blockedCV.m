function [X_train, Y_train, X_test, Y_test, cv_blocked] = ...
    ml_blockedCV(X,Y, Blocks,Stratify, Kfolds)
    [foldIdx, cruiseFold] = makeCruiseMonthFolds(Blocks, Stratify, Kfolds);

    %Test and training
    X_train = X(:, :);
    Y_train = Y(:);
    X_test = [];
    Y_test = [];

    cv_blocked = cvpartition("CustomPartition", foldIdx);
    nFolds = cv_blocked.NumTestSets;
end