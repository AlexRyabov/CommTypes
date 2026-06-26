%Feature selection
clearvars
proj_Initialize;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
%LM_Model = 'Boost_all'; %boosted regression with all predictors
%LM_Model = 'Boost_noimp'; %boosted regression no imputed predict-
%LM_Model = 'Boost_loc'; %use date and location
%LM_Model = 'Boost_env'; %use date and location
%LM_Model = 'Boost_env_sat'; %CTD, nutrients, satellite
%LM_Model = 'Boost_regular'; %regularly sampled data and most important (based on boosted )
%LM_Model = 'Boost_regular_no_depthmonth'; %regularly sampled data and most important, without depth and month info
%LM_Model = 'Boost_ctd_sat'; %regularly sampled ctd and satellite data
%LM_Model = 'Boost_ctd_sat_no_depth,doy' %regularly sampled ctd and satellite data, no depth and doy
%LM_Model = 'Boost_sat_depth'; %satellite data and depth
%LM_Model = 'Boost_selected'; %start from regular sampled and then remove less important
%LM_Model = 'Boost_most'; %satellite, imputed nutrients, depth
%LM_Model = 'Boost_best_varselect';  %after variable selection
%LM_Model = 'Boost_best_varselect+depth'  %after variable selection + depth
%LM_Model = 'Temp,Nutr,depth,Doy'  %Temp, nutr, depth and doy
%LM_Model = 'Boost_RegAndImputed'; LM_Model_title = 'Env, depth, doy'; %regularly sampled ice data, imputed environment, depth and doy
LM_Model = 'Boost_RegAndImputedNew'; LM_Model_title = 'Env, depth, doy'; %regularly sampled ice data, new imputed environment, depth and doy

%LM_Model = 'Boost_RegAndImputedNoDoy'; LM_Model_title = 'Env, depth'; %regularly sampled ice data, imputed environment, depth no doy

Fram_0400_Init;

%Class model ID
ClassifID = 'ECOC_SVM';
%ClassifID = 'RF';

trainParams = {};
classes = categories(Y_train);
K = numel(classes);

RemoveDoY = false;
RemoveDepth = false;

if RemoveDoY
    ind = strcmp('Month_doy', predictors);
    if any(ind)
        predictors = predictors(~ind);
        X_train = X_train(:, predictors);
        if ~isempty(X_test)
            X_test = X_test(:, predictors);
        end
        X = X(:, predictors);
    end
end

if RemoveDepth
    ind = strcmp('Depth', predictors);
    if any(ind)
        predictors = predictors(~ind);
        X_train = X_train(:, predictors);
        if ~isempty(X_test)
            X_test = X_test(:, predictors);
        end
        X = X(:, predictors);
    end
end

switch ClassifID
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = [1, 2, 3];
        HP_grid.boxconstr  = [100, 150, 200];

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
predictors(mask_seq2)


%% Refit hyper parameters
switch ClassifID
    case 'ECOC_SVM'
        %% hyperparameter grid
        HP_grid.sigmas     = logspace(0, 1, 10);
        HP_grid.boxconstr  = logspace(-1, 2, 12);
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
switch ClassifID
    case 'ECOC_SVM'
        Mdl1 = ml_train_func(X_train(:, mask_seq2), Y_train, bestParams,  'FitPosterior', true);
    otherwise
        Mdl1 = ml_train_func(X_train(:, mask_seq2), Y_train, bestParams);
end

[cvErr, CVMdl] = mlTrainAndCVLoss(X_train(:, mask_seq2), Y_train, cv_inner, ml_train_func, bestParams);
disp(['CV err ' NS(cvErr)])

% Suppose you trained cross-validated model, e.g.:
% CVMdl = fitcecoc(X, Y, 'Learners', tSVM, 'CVPartition', cv);

[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,400 400], 'normal', 'on');

% 1) Out-of-fold predictions
yHat = kfoldPredict(CVMdl);   % each point predicted by a model that did not see it

% 2) Confusion matrix (counts)
C = confusionmat(Y_train, yHat);

% 3) Row-normalized confusion matrix (per-class mixing)
Cnorm = C ./ sum(C, 2);

% 4) Nice plot
%confusionchart(Y_train, yHat);      % by default shows counts
confusionchart(Y_train, yHat, 'Normalization','row-normalized');  % mixing per true class
title(LM_Model_title)
%Save the algorithm

% save(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']), 'Mdl1', 'ClassifID', ...
%     'mask_seq2', 'cv_inner', 'bestParams', 'cvErr', 'LM_Model_title');

%Fram_0435_EnvType_PartDepend

%%
