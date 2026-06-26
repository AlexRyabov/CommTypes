%train



%% hyperparameter tuning
%% DATA

Fram_040_Init


%% GRID
% Restrict to the inner-CV training subset defined by you

% Fix class order across folds
classes = categories(Y_train);

if 1 %make parameter tuning

    %% YOUR CUSTOM BLOCK PARTITION (provided logic)
    nFolds = cv_inner.NumTestSets;

    %%HYPERPARAMETER GRID
    grid.sigmas = logspace(log10(4), log10(8), 10);
    grid.boxconstr  = logspace(1, 2.5, 10);

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

    %%
    errMat = NaN(height(tParams), 1);
    parfor iRow = 1:height(tParams)
        iRow
        Sigma = tParams.sigmas(iRow);
        C = tParams.boxconstr(iRow);
        t = templateSVM( ...
            'KernelFunction', 'gaussian', ...  % RBF
            'KernelScale',    Sigma, ...         % sigma
            'BoxConstraint',  C, ...           % C
            'Standardize',    true, ...        % 
            'Type',           'classification');

        CVMdl = fitcecoc( ...
            X_train, Y_train, ...
            'Learners',     t, ...
            'CVPartition',  cv_inner, 'Coding','onevsall');
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
    

%%
    %%/ best =
    %
    %   struct with fields:
    %
    %         err: 0.5749
    %      NTrees: 100
    %     boxconstr: 10
    %        Mtry: 5
    %
    %
    fprintf('Best (custom block-CV): Err=%.4f, sigmas=%d, boxconstr=%d\n', ...
        best.err, best.sigmas, best.boxconstr);
end
%% FINAL MODEL on the entire inner training subset (optional)

% best = 
%           err: 0.4160
%     boxconstr: 14.6780
%        sigmas: 5.0397


       % t = templateSVM( ...
       %      'KernelFunction', 'gaussian', ...  % RBF
       %      'KernelScale',    best.sigmas, ...         % sigma
       %      'BoxConstraint',  best.boxconstr, ...           % C
       %      'Standardize',    true, ...        % 
       %      'Type',           'classification');


       nTrees = 200;
minLeaf = 10;
Mtry = 3;

 cvError = blockCVError(nTrees,  minLeaf, Mtry)
 disp(["CVError " + cvError])

Mdl1 = TreeBagger(nTrees, X_train, Y_train, ...
    'Method','classification', ...
    'MinLeafSize', minLeaf, ...
    'NumPredictorsToSample', Mtry, ...
    'Surrogate','off', ...
    'OOBPrediction','on', 'OOBPredictorImportance', 'on');
% 
t = templateEnsemble('Bag', 100, 'Tree', ...
    'Type', 'Classification', ...
    'Learners', templateTree('MaxNumSplits', 50, 'MinLeafSize', minLeaf));



 


% OOB permutation importance (Δ OOB classification error per feature)


Fram_0490_PlotResultsBag

