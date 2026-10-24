function err = mlTrainAndLoss(X, Y, X_Test, Y_Test, funcTrain, params, varargin)
%Train the model and find its loss 
mld = funcTrain(X, Y, params, varargin{:});
err= loss(mld, X_Test, Y_Test);
end
