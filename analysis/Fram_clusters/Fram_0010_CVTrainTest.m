%Generate cross validation, test and train table 


if HGP.ML.CVBlocked
   [foldIdx, cruiseFold] = makeCruiseMonthFolds(tStPredictors.Cruise, tStPredictors.Month, HGP.ML.KFolds);
   cv_inner = cvpartition("CustomPartition", foldIdx);
else 
   cv_inner = cvpartition(MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
end

%Test and training
X_train = X(:, :);
Y_train = Y(:);
X_test = [];
Y_test = [];