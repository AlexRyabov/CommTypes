%train



%% hyperparameter tuning
%% DATA

Fram_0400_Init


%% GRID
% Restrict to the inner-CV training subset defined by you

% Fix class order across folds
classes = categories(Y_train);

ix = @(c) find(strcmp(classes, string(c)));

K = numel(classes);
Cost = ones(K) - eye(K);                    % base cost
Cost(ix(1), [ix(3)]) = 2;       % penalties for leaking 1
Cost(ix(6), [ix(3) ix(9)]) = 2;       % penalties for leaking  6
Cost(ix(8), [ix(2) ix(9)]) = 2;       % penalties for leaking 8

others = setdiff(1:K, ix(9));
Cost(others, ix(9)) = 2;   % tune 2.0–4.0 on validation

%% make parameter tuning
if 1 
   nFolds = cv_inner.NumTestSets;

    %%HYPERPARAMETER GRID
    grid.sigmas = logspace(log10(4), log10(20), 10);
    grid.boxconstr  = logspace(0.1, 2.5, 10);

    P = numel(predictors);

    % %% EVALUATION FUNCTION (misclassification rate averaged over custom folds)
    % blockCVError = @(numTrees, boxconstr, mtry) mean( ...
    %     arrayfun(@(k) foldErr(k, numTrees, boxconstr, mtry, X_train, Y_train, cv_inner, classes), 1:nFolds) );

    %%TUNING + PLOTTING
    numSgm = numel(grid.sigmas);
    numC = numel(grid.boxconstr);

    ErrMaps = nan(numSgm, numC);
    best.err = inf;

    tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');

    tParams = table();
    iRow = 1;
    for iSgm = 1:numSgm
        for iC = 1:numC
            tParams.sigmas(iRow) = grid.sigmas(iSgm);
            tParams.boxconstr(iRow) = grid.boxconstr(iC);
            tParams.iSgm(iRow) = iSgm;
            tParams.iC(iRow) = iC;
            iRow = iRow +1;
        end
    end

    %
    errMat = NaN(height(tParams), 1);
    parfor iRow = 1:height(tParams)
        iRow
        Sigma = tParams.sigmas(iRow);
        C = tParams.boxconstr(iRow);


        tSVM = templateSVM('KernelFunction','rbf','Standardize',true,  'KernelScale',   ...
            Sigma,  'BoxConstraint',  C,  'Type',  'classification');
        CVMdl = fitcecoc(X_train, Y_train, ...
            'Learners', tSVM, ...
            'Coding', 'onevsone', ...
            'Cost', Cost, ...
            'CVPartition',  cv_inner);


        errMat(iRow) = kfoldLoss(CVMdl);
    end

    tParams.errMat = errMat;
    %plot the results
    for iRow = 1:height(tParams)
        ErrMaps(tParams.iSgm(iRow), tParams.iC(iRow)) = tParams.errMat(iRow);
    end

    figure(1)
    clf
    nexttile;
    [BC, SGM] = meshgrid(grid.boxconstr, grid.sigmas);
    contourf(BC, SGM, ErrMaps, 'ShowText','on'); colorbar;
    f_Set_scales({'logx', 'logy'})
    f_Lbls('C', 'sigma')
    [~, indMin] = min( tParams.errMat);
    hold on
    scatter(tParams.boxconstr(indMin), tParams.sigmas(indMin), 30, 'filled')
    hold off


    best.boxconstr = tParams.boxconstr(indMin);
    best.sigmas = tParams.sigmas(indMin);
    best.err = tParams.errMat(indMin);



    %
    fprintf('Best (custom block-CV): Err=%.4f, sigmas=%d, boxconstr=%d\n', ...
        best.err, best.sigmas, best.boxconstr);
end
%% FINAL MODEL on the entire inner training subset (optional)

% 
% best.err= 0.3141;
% best.boxconstr=  1.2589;
% best.sigmas=4;


t = templateSVM( ...
    'KernelFunction', 'rbf', ...  % RBF
    'BoxConstraint',  best.boxconstr, ...           % C
    'KernelScale',  best.sigmas, ...           % C
    'Standardize',    true, ...        %
    'Type',           'classification');

Mdl1 = fitcecoc( ...
    X_train, Y_train, ...
    'Learners',     t, 'Coding','onevsall',  'Cost', Cost,  'FitPosterior', true);

%%  Add penalties for class 6 and 8


tSVM = templateSVM('KernelFunction','rbf','Standardize',true, 'KernelScale', 'auto');
Mdl1 = fitcecoc(X_train, Y_train, ...
    'Learners', tSVM, ...
    'Coding', 'onevsone', ...
    'Cost', Cost,  ...
    'OptimizeHyperparameters', {'BoxConstraint'}, ...
    'HyperparameterOptimizationOptions', struct( 'CVPartition',  cv_inner,'MaxObjectiveEvaluations',40,'UseParallel',true));


CV = crossval(Mdl1,'CVPartition',  cv_inner);
trLoss = resubLoss(Mdl1)
cvLoss = kfoldLoss(CV)

%% linear SVM 
%%  Add penalties for class 6 and 8


tSVM = templateSVM('KernelFunction','linear','Standardize',true);
Mdl1 = fitcecoc(X_train, Y_train, ...
    'Learners', tSVM, ...
    'Coding', 'onevsone', ...
    'Cost', Cost,  ...
    'OptimizeHyperparameters', {'BoxConstraint'}, ...
    'HyperparameterOptimizationOptions', struct( 'CVPartition',  cv_inner,'MaxObjectiveEvaluations',40,'UseParallel',true));


CV = crossval(Mdl1,'CVPartition',  cv_inner);
trLoss = resubLoss(Mdl1)
cvLoss = kfoldLoss(CV)


    %  BoxConstraint
    % _____________
    % 
    %     0.27236   



%%
% % Rank features quickly (filter stage)
% [idx,score] = fscmrmr(X_train, Y_train);
% k0 = min(50, numel(idx));                % keep top-k0 as a prefilter
% X_train0 = X_train(:, idx(1:k0));

%% 
% === Outer CV for feature selection ===
% opts = statset('UseParallel', true);     % speed up
% BoxConstraint =  0.3015;
% % Criterion: train ECOC+SVM with inner CV (tunes C) and score on hold-out
% critfun = @(Xtr, Ytr, Xte, Yte) localLossECOC(Xtr, Ytr, Xte, Yte, Cost, BoxConstraint);
% 
% % Run SFS on the prefiltered feature set (optional but speeds up a lot)
% [inModel, history] = sequentialfs(critfun, X_train{:, :}, Y_train{:}, ...
%     'CV', cv_inner, 'direction', 'forward', 'nfeatures', 20, 'options', opts);
% 
% % Selected features in original X_train coordinates:
% selectedIdx = idx(inModel);              % map back from prefiltered indices
% X_train_sel = X_train(:, selectedIdx);
% 
% % Train final model on selected features
% tSVM = templateSVM('KernelFunction','rbf','Standardize',true,'KernelScale','auto');
% Mdl = fitcecoc(X_train_sel, Y_train, 'Learners', tSVM, 'Coding','onevsone','Cost',Cost);
% 


%% classification with RUSBoost
%
% tTree = templateTree('MaxNumSplits',1,'MinLeafSize',5);
% tRUS  = templateEnsemble('RUSBoost',300,tTree,'LearnRate',0.05);
% Mdl1 = fitcecoc(X_train, Y_train, 'Learners', tRUS, 'Coding','onevsone', 'Cost', Cost);

%% plot results


Fram_0490_PlotResultsBag

