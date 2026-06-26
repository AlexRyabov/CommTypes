function loss = ml_cvLossNN(params, X, Y, cv_blocked)
% params: table from bayesopt with field:
%   NumHidden   - number of hidden units (1..20)
% X, Y: data
% cv_blocked: cvpartition object

    % Cross-validated regression network with 1 hidden layer
  MdlCV = fitrnet( ...
        X, Y, ...
        'LayerSizes', params.NumHidden, ...   % one hidden layer
        'Lambda',     params.Lambda, ...      % L2 regularization
        'Standardize', true, ...
        'CVPartition', cv_blocked);
  
    % k-fold loss (MSE by default)
    loss = kfoldLoss(MdlCV);
end
