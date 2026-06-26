%selection of variables in random forest 



%% hyperparameter tuning
%% DATA
% T is your table.
% Assumes:
%   - response column: 'Y' (class labels)
%   - blocking column: 'Year' (one block = one year)
%   - all other columns are predictors

Fram_0400_Init

%% GRID
% Restrict to the inner-CV training subset defined by you

% Fix class order across folds
classes = categories(Y_train);
NumTrees = 300;
MinLeaf = 10;
nFolds = cv_inner.NumTestSets;
if 1 %make variable selection 
    min_n_predictors = 2; %minimal number of predictors 
    max_n_predictors = length(predictors); %max number of predictors 
    tResults = table;
    tResults.predictors = cell(max_n_predictors, 1);
    tResults.cv_error = NaN(max_n_predictors, 1);
    cur_predictors = predictors;
    %%
    for iP = 1:length(predictors)-min_n_predictors 
        %fit with cur_predictors
        tResults.predictors{iP} = cur_predictors;
        P = length(cur_predictors);
         tResults.P(iP) = P;
        tResults.Mtry(iP) =  3; %max(1, ceil(P/3));
        X_train_cur =  X_train(:, cur_predictors);
       %%EVALUATION FUNCTION (misclassification rate averaged over custom folds)
        blockCVError = @(numTrees, minLeaf, mtry) mean( ...
            arrayfun(@(k) foldErr(k, numTrees, minLeaf, mtry,X_train_cur, Y_train, cv_inner, classes), 1:nFolds) );
        tResults.cv_error(iP) = blockCVError( NumTrees,  MinLeaf,  tResults.Mtry(iP) );

        Mdl1 = TreeBagger(NumTrees, X_train_cur, Y_train, ...
            'Method','classification', ...
            'MinLeafSize', MinLeaf, ...
            'NumPredictorsToSample', tResults.Mtry(iP), ...
            'Surrogate','on', ...
            'OOBPrediction','on', 'OOBPredictorImportance', 'on');
        pred_importance = Mdl1.OOBPermutedPredictorDeltaError;
        [~, idx] = sort(pred_importance, "descend");
        figure(1)
        nexttile(1)
        bar(cleanUndersc(X_train_cur.Properties.VariableNames(idx)),(pred_importance(idx)));
        nexttile(2)
        plot(tResults.P, tResults.cv_error)
        drawnow;
        new_pred = ~(idx(end) == 1:P);
        cur_predictors = cur_predictors(new_pred);
        tResults
    end
 
end

[~, ind] = min(tResults.cv_error + tResults.P/100);
%[~, ind] = min(tResults.cv_error );
predictors = tResults.predictors{ind};
X_train =  X_train(:, predictors);
X_test =  X_test(:, predictors);
Mdl1 = TreeBagger(100,  X(:, predictors), Y, ...
    'Method','classification', ...
    'MinLeafSize', MinLeaf, ...
    'NumPredictorsToSample',  tResults.Mtry(ind), ...
    'Surrogate','on', ...
    'OOBPrediction','on', 'OOBPredictorImportance', 'on');
save(fullfile(HGP.Data.Path2ResultsData, ['trainedRF_' LM_Model '.mat']), 'Mdl1');

% OOB permutation importance (Δ OOB classification error per feature)

Fram049_PlotResultsBag

%bestPred = 

%% ---- Local function (place at end of file) ----
function e = foldErr(k, numTrees, minLeaf, mtry, Xtr, Ytr, cv_inner, classes)
% Binary masks provided by your custom cvpartition
trMask = cv_inner.training(k);
teMask = cv_inner.test(k);

% Train Random Forest on training block
M = TreeBagger(numTrees, Xtr(trMask,:), Ytr(trMask), ...
    'Method','classification', ...
    'MinLeafSize', minLeaf, ...
    'NumPredictorsToSample', mtry, ...
    'Surrogate','off', ...
    'OOBPrediction','off');

% Predict on held-out block
Yhat = predict(M, Xtr(teMask,:));     % cellstr
Yhat = categorical(Yhat, classes);    % align class order
e = mean(Yhat ~= Ytr(teMask));        % misclassification rate
end