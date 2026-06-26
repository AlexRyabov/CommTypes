%Load data and define parameters of LM model to predict the class based on 
%Lon, Month and Depth
%or Temperature, Month and Depth 

proj_Initialize;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));

%% select the data for temperature, depth and location
predictors_loc = { 'temp_imputed', 'Month_doy', 'Depth'};

predictors = predictors_loc;
%LM_Model = 'ECOC_Location'; %location
LM_Model = 'ECOC_Location_randCV'; %temperature, stratification by month, no blocks by year

X = tStPredictors(:, predictors);

Y = categorical(tDMC.ClustID);  %select nonnan Y

[fg, MaxFigure] = addFindFigure('Present data', MaxFigure, [10, 20,400 400], 'normal', 'on');
t = nexttile
ind_missed_data = ismissing(X);
notMissed = sum(~ind_missed_data);
bar(predictors, notMissed)
%predictors(notMissed >700)';
ind_missed_rows  = sum(ind_missed_data, 2) >0;
X = X(~ind_missed_rows, :);
Y = Y(~ind_missed_rows);
if height(X) < 600
    error('too many rows were removeded')
end

% for blocked cross validation
Cruise = tStPredictors.Cruise(~ind_missed_rows);
MonthNum = tStPredictors.MonthNum(~ind_missed_rows);




if HGP.ML.CVBlocked
    [foldIdx, cruiseFold] = makeCruiseMonthFolds(Cruise, MonthNum, HGP.ML.KFolds);
    cv_inner = cvpartition("CustomPartition", foldIdx);
else  %location, stratification by month, no blocks by year
    cv_inner = cvpartition(MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
end
nFolds = cv_inner.NumTestSets;

%Test and training
X_train = X(:, :);
Y_train = Y(:);
X_test = [];
Y_test = [];





