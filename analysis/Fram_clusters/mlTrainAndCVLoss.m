function [cvErr, CVMdl] = mlTrainAndCVLoss(X, Y, cv, funcTrain, params, varargin)

CVMdl = funcTrain(X, Y, params, 'CVPartition', cv, varargin{:});

cvErr = kfoldLoss(CVMdl);
end
