%Feature selection

Fram_0310_MLLocation_Init;

%Class model ID
ClassifID = 'ECOC_SVM';
%ClassifID = 'RF';

trainParams = {};
classes = categories(Y_train);

%%

switch ClassifID
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = logspace(log10(0.1), log10(30), 10);
        HP_grid.boxconstr  = logspace(-1, 2, 11);

        ix = @(c) find(strcmp(classes, string(c)));
        %
        K = numel(classes);
        GroupCounts = grpstats(table(Y_train, 'VariableNames', {'Y'}), 'Y', [ 'numel']);
        GroupCounts.Weight = 1./GroupCounts.GroupCount;
        % [~, loc] = ismember(Y_train, GroupCounts.Y);
        % obsWeights = GroupCounts.Weight(loc);

    
        Cost = ones(K) - eye(K);          % base 0/1 cost
        % define the cost 
        for k = 1:K
            Cost(k, :) = Cost(k, :) * GroupCounts.Weight(k);
        end

        ECOC_Options = {'Cost', Cost};
        %ECOC_Options = {'Weights', obsWeights};

        %handle of the training function ml_train_func(X, Y, params, Name, Value....)
        %ml_train_func = @(X, Y, params, varargin) ml_train_ECOC(X, Y, params, ECOC_Options{:}, varargin{:});
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

        % trainParams = {'Surrogate','on',  'OOBPrediction','on'};

        %handle of the training function ml_train_func(X, Y, params, Name, Value....)
        ml_train_func = @ml_train_RF;

        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train, Y_train, cv_inner, ml_train_func, HP_grid, trainParams{:});
        bestParams
end

%% Get test loss for sequential fit
fsCritFun = @(X_tr, Y_tr,X_te, Y_te) ...
    mlTrainAndLoss(X_tr, Y_tr, X_te, Y_te, ml_train_func, bestParams);


[mask_seq, history] = sequentialfs(fsCritFun, X_train{:, :}, Y_train, 'CV', cv_inner,"Direction", "backward", "Options",statset("Display","iter", "UseParallel",true));
predictors(mask_seq)


%% Refit hyper parameters
switch ClassifID
    case 'ECOC_SVM'
        %% hyperparameter grid

        HP_grid.sigmas     = logspace(0, 1, 10);
        HP_grid.boxconstr  = logspace(-1, 2, 12);
        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train(:, mask_seq), Y_train, cv_inner,  ml_train_func, HP_grid);

        figure(fg);
        nexttile(2)
        Sigmas   = paramMesh{1};   % MaxNumSplits
        BC  = paramMesh{2};   % NumLearningCycles

        contourf(HP_grid.boxconstr, HP_grid.sigmas, CVloss, 'ShowText', 'on'); colorbar;
        f_Set_scales({'logx', 'logy'})
        xlabel('BC')
        ylabel('Sigma')
        hold on
        scatter(bestParams.boxconstr, bestParams.sigmas, 30, 'red', 'filled');
        hold off
        title('Best params, selected features')
        drawnow

    case 'RF'
        %% hyperparameter grid
        HP_grid.NumTrees = [50, 100, 200, 300];
        HP_grid.MinLeafSize  = [5:5:20];
        P = numel(predictors);
        HP_grid.NumPredictorsToSample = unique([max(1, round(sqrt(P))), ceil(P/3), P]); % extern

        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train(:, mask_seq), Y_train, cv_inner, ml_train_func, HP_grid, trainParams{:});


    otherwise
        error('ClassifID not defined')
end


%% Try to select features from the selected set
%Get test loss for sequential fit
fsCritFun = @(X_tr, Y_tr,X_te, Y_te) ...
    mlTrainAndLoss(X_tr, Y_tr, X_te, Y_te, ml_train_func, bestParams);
selectedPredictors = predictors(mask_seq);
[mask_seq_tmp, history] = sequentialfs(fsCritFun, X_train{:, mask_seq}, Y_train, 'Direction', 'forward', 'CV', cv_inner, "Options",statset("Display","iter"));
selectedPredictors2 = selectedPredictors(mask_seq_tmp)
mask_seq2 = ismember(predictors, selectedPredictors2);


%% Refit hyper parameters
switch ClassifID
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = logspace(0, 1, 10);
        HP_grid.boxconstr  = logspace(-1, 2.5, 12);
        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train(:, mask_seq2), Y_train, cv_inner,  ml_train_func, HP_grid);

        figure(fg);
        nexttile(3)
        Sigmas   = paramMesh{1};   % MaxNumSplits
        BC  = paramMesh{2};   % NumLearningCycles

        contourf(HP_grid.boxconstr, HP_grid.sigmas, CVloss, 'ShowText', 'on'); colorbar;
        f_Set_scales({'logx', 'logy'})
        xlabel('BC')
        ylabel('Sigma')
        hold on
        scatter(bestParams.boxconstr, bestParams.sigmas, 30, 'red', 'filled');
        hold off
        title('Best params, selected features')
        drawnow

    case 'RF'
        %% hyperparameter grid
        HP_grid.NumTrees = [25, 50, 100];
        HP_grid.MinLeafSize  = [2 5, 10, 15, 20, 25];
        P = numel(predictors);
        HP_grid.NumPredictorsToSample = unique([max(1, round(sqrt(P))), ceil(P/3), P]); % extern

        [bestParams, CVloss, paramMesh] = gridSearchCV(X_train(:, mask_seq2), Y_train, cv_inner, ml_train_func, HP_grid, trainParams{:});
        bestParams

    otherwise
        error('ClassifID not defined')
end

%%

Mdl1 = ml_train_func(X_train(:, mask_seq2), Y_train, bestParams,  'FitPosterior', true);

[cvErr, CVMdl] = mlTrainAndCVLoss(X_train(:, mask_seq2), Y_train, cv_inner, ml_train_func, bestParams);
disp(['CV err ' NS(cvErr)])

% Suppose you trained cross-validated model, e.g.:
% CVMdl = fitcecoc(X, Y, 'Learners', tSVM, 'CVPartition', cv);

[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');

% 1) Out-of-fold predictions
yHat = kfoldPredict(CVMdl);   % each point predicted by a model that did not see it

% 2) Confusion matrix (counts)
C = confusionmat(Y_train, yHat);

% 3) Row-normalized confusion matrix (per-class mixing)
Cnorm = C ./ sum(C, 2);

% 4) Nice plot
%confusionchart(Y_train, yHat);      % by default shows counts
confusionchart(Y_train, yHat, 'Normalization','row-normalized');  % mixing per true class

%Save the algorithm 
save(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']), 'Mdl1');