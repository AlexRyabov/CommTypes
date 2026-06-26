%train random forest classification



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

if 0 %make parameter tuning

    %% YOUR CUSTOM BLOCK PARTITION (provided logic)
    nFolds = cv_inner.NumTestSets;

    %% HYPERPARAMETER GRID
    grid.NumTrees = [50:50:500];
    grid.MinLeaf  = [5:1:20];

    P = numel(predictors);
    grid.Mtry = unique([max(1, round(sqrt(P))), ceil(P/3), P]); % external loop
    %grid.Mtry = 3:10; % external loop

    %% EVALUATION FUNCTION (misclassification rate averaged over custom folds)
    blockCVError = @(numTrees, minLeaf, mtry) mean( ...
        arrayfun(@(k) foldErr(k, numTrees, minLeaf, mtry, X_train, Y_train, cv_inner, classes), 1:nFolds) );

    %% TUNING + PLOTTING
    numN = numel(grid.NumTrees);
    numL = numel(grid.MinLeaf);
    numV = numel(grid.Mtry);

    ErrMaps = nan(numL, numN, numV);  % [MinLeaf x NumTrees x Mtry]
    best.err = inf;

    tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');

    tParams = table();
    iRow = 1;
    for vIdx = 1:numV
        mtry = grid.Mtry(vIdx);  % external variable as requested

        % Compute error surface (MinLeaf × NumTrees) for this Mtry
        errMat = nan(numL, numN);
        for iL = 1:numL
            for iN = 1:numN
                tParams.mtry(iRow) = mtry;
                tParams.NumTrees(iRow) = grid.NumTrees(iN);
                tParams.MinLeaf(iRow) = grid.MinLeaf(iL);
                iRow = iRow +1;
            end
        end
    end

    errMat = NaN(height(tParams), 1);
    parfor iRow = 1:height(tParams)
        errMat(iRow) = blockCVError( tParams.NumTrees(iRow),  tParams.MinLeaf(iRow),  tParams.mtry(iRow));
    end
    tParams.errMat = errMat;
    %%
    %plot the results
    for vIdx = 1:numV
        mtry = grid.Mtry(vIdx)  % external variable as requested
        % Compute error surface (MinLeaf × NumTrees) for this Mtry
        errMat = nan(numL, numN);
        for iL = 1:numL
            for iN = 1:numN
                ind = tParams.mtry == mtry & tParams.NumTrees == grid.NumTrees(iN) & tParams.MinLeaf== grid.MinLeaf(iL)
                errMat(iL, iN) = tParams.errMat(ind);
            end
        end
        ErrMaps(:,:,vIdx) = errMat;
    end
    minVal= min(ErrMaps(:));
    maxVal= max(ErrMaps(:));
    %%
    figure(1)
    clf
    best.err = Inf;
    for vIdx = 1:numV
        mtry = grid.Mtry(vIdx)  % external variable as requested
        errMat = ErrMaps(:,:,vIdx);
        % Best point on this surface
        [minVal, linIdx] = min(errMat(:));
        [iLbest, iNbest] = ind2sub(size(errMat), linIdx);
        Nbest = grid.NumTrees(iNbest);
        Lbest = grid.MinLeaf(iLbest);

        % Track global best
        if minVal < best.err
            best.err = minVal;
            best.NTrees  = Nbest;
            best.MinLeaf = Lbest;
            best.Mtry    = mtry;
        end

        % ----- Contour plot for this Mtry -----
        nexttile;
        [NN, LL] = meshgrid(grid.NumTrees, grid.MinLeaf);
        contourf(NN, LL, errMat, 'ShowText','on'); colorbar;
        hold on;
        caxis([minVal, maxVal]);
        plot(Nbest, Lbest, 'wo', 'MarkerFaceColor','k', 'MarkerSize',6);
        hold off;
        xlabel('NumTrees'); ylabel('MinLeafSize');
        title(sprintf('Block-CV error (custom folds), Mtry = %d', mtry));
        set(gca,'YDir','normal');
        drawnow
    end
    %%/ best =
    %
    %   struct with fields:
    %
    %         err: 0.5749
    %      NTrees: 100
    %     MinLeaf: 10
    %        Mtry: 5
    %
    %
    fprintf('Best (custom block-CV): Err=%.4f, NumTrees=%d, MinLeaf=%d, Mtry=%d\n', ...
        best.err, best.NTrees, best.MinLeaf, best.Mtry);
end
%% FINAL MODEL on the entire inner training subset (optional)

    nFolds = cv_inner.NumTestSets;
    %% EVALUATION FUNCTION (misclassification rate averaged over custom folds)
    blockCVError = @(numTrees, minLeaf, mtry) mean( ...
        arrayfun(@(k) foldErr(k, numTrees, minLeaf, mtry, X, Y, cv_inner, classes), 1:nFolds) );


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


% t = templateTree('MinLeafSize',minLeaf, ...
%                  'NumVariablesToSample', Mtry, ...
%                  'Surrogate','off');
% 
% Mdl1 = fitcensemble(X_train, Y_train, ...
%     'Method','Bag', ...
%     'Learners', t, ...
%     'NumLearningCycles', nTrees);

% OOB permutation importance (Δ OOB classification error per feature)

Fram_0490_PlotResultsBag

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