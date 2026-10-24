function [best, CVloss, paramMesh] = ...
    gridSearchCV(X, Y, cv, funcTrain, paramGrid, varargin)
%GRIDSEARCHCV Hyperparameter grid search with cross-validation.
%
%   [BEST, CVLOSS, PARAMMESH] = GRIDSEARCHCV(X, Y, CV, PARAMGRID, FITFCN)
%   performs an exhaustive grid search over hyperparameters defined in
%   PARAMGRID. For each combination, FITFCN is called as:
%
%       loss = FITFCN(X, Y, CV, params, ....)
%
%   where PARAMS is a struct with one field per hyperparameter.
%
%   INPUT:
%       X         - Predictor matrix (observations x features).
%       Y         - Response vector (class labels or regression targets).
%       cv        - cvpartition object defining cross-validation folds.
%       paramGrid - Struct with one field per hyperparameter, for example:
%                   paramGrid.MaxNumSplits      = [5 10 20 50];
%                   paramGrid.NumLearningCycles = [50 100 200];
%       fitFcn    - Function handle that returns cross-validated loss
%                   for a given parameter struct:
%                   loss = fitFcn(X, Y, cv, params, ...)
%
%   OUTPUT:
%       best      - Struct with the best hyperparameter values and field:
%                   best.err = minimal cross-validated loss.
%       CVloss    - N-D array of losses with the same size as the
%                   hyperparameter grids in PARAMMESH.
%       paramMesh - Cell array of N-D grids produced by NDGRID for each
%                   hyperparameter in PARAMGRID.
%
%   EXAMPLE (fitcensemble with AdaBoost and tree learners):
%
%       function err = localLossEnsemble(X, Y, cv_inner, params, method)
%
%           % params.MaxNumSplits, params.NumLearningCycles
%           tTree = templateTree( ...
%               'MaxNumSplits', params.MaxNumSplits);
%
%           CVMdl = fitcensemble(X, Y, ...
%               'Method',           method, ...
%               'NumLearningCycles', params.NumLearningCycles, ...
%               'Learners',         tTree, ...
%               'CVPartition',      cv_inner);
%
%           err = kfoldLoss(CVMdl);
%       end
%
%       % Hyperparameter grid
%       grid.MaxNumSplits      = [5 10 20 50];
%       grid.NumLearningCycles = [50 100 200];
%
%       method = 'AdaBoostM1';
%
%       fitFcnEnsemble = @(X, Y, cv_inner, params) ...
%           localLossEnsemble(X, Y, cv_inner, params, method);
%
%       [best, CVloss, paramMesh] = ...
%           gridSearchCV(X_train, Y_train, cv_inner, grid, fitFcnEnsemble);
%
%       % Example plotting for a 2D grid:
%       % MS   = paramMesh{1};   % MaxNumSplits
%       % NLC  = paramMesh{2};   % NumLearningCycles
%       % contourf(NLC, MS, CVloss, 'ShowText', 'on'); colorbar;
%       % xlabel('NumLearningCycles'); ylabel('MaxNumSplits');


% Extract hyperparameter names and number of dimensions
paramNames = fieldnames(paramGrid);
nParams    = numel(paramNames);

% Create full factorial N-D grid of all hyperparameter combinations
gridVectors = struct2cell(paramGrid);          % {grid.MaxNumSplits, grid.NumLearningCycles, ...}
[paramMesh{1:nParams}] = ndgrid(gridVectors{:});

% Total number of combinations (all paramMesh{i} share the same size)
nComb = numel(paramMesh{1});

% Build table of all combinations in vector form
paramList = table();
for i = 1:nParams
    vals = paramMesh{i}(:);                    % flatten to column vector
    paramList.(paramNames{i}) = vals;
end

% Preallocate loss vector
errMat = NaN(nComb, 1);

parfor i = 1:nComb
    % Build params struct for the current combination
    [i, nComb]
    params = struct();
    for j = 1:nParams
        pname = paramNames{j};
        params.(pname) = paramList.(pname)(i);
    end

    % FITFCN must return cross-validated loss for the given params
    CVMdl = funcTrain(X, Y, params, 'CVPartition', cv, varargin{:});
    errMat(i) = kfoldLoss(CVMdl);

end

paramList.errMat = errMat;

% Reshape vector of losses back to the N-D grid shape
gridSize = size(paramMesh{1});    % e.g. [numMaxNumSplits, numNumLearningCycles]
CVloss   = reshape(errMat, gridSize);
% Now CVloss(i1, i2, ..., iN) corresponds to
% paramMesh{1}(i1,...,iN), paramMesh{2}(i1,...,iN), ..., paramMesh{N}(i1,...,iN).

% Find best combination in vector form
[bestErr, indMin] = min(errMat);

best     = struct();
best.err = bestErr;
for j = 1:nParams
    pname        = paramNames{j};
    best.(pname) = paramList.(pname)(indMin);
end
end
