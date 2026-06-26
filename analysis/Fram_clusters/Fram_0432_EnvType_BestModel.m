%Feature selection
clearvars
proj_Initialize;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
LM_Model = 'Boost_best_varselect'; LM_Model_title = 'temp_imputed, depth, doy'; %Best predictors Temp, Depth and doy 


Fram_0400_Init;

%Class model ID
ClassifID = 'ECOC_SVM';

trainParams = {};
classes = categories(Y_train);
K = numel(classes);

mask_seq2 = true(size(predictors));

%% fit rbf SVM model 
%template
tSVM = templateSVM( ...
    'KernelFunction', 'rbf', ...
    'Standardize',    true, ...
    'Type',           'classification');
%options for automatic fitting of paramters 
opt = struct( ...
    'CVPartition', cv_inner, ...
    'Optimizer', 'bayesopt', ...
    'AcquisitionFunctionName','expected-improvement-plus', ...
    'ShowPlots', true, ...
    'Verbose', 1, ...
    'MaxObjectiveEvaluations', 30);

%parameter optimization 
[Mdl] = fitcecoc(X_train, Y_train, ...
    'Learners', tSVM, ...
    'Coding', 'onevsone', ...
    'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ...
    'HyperparameterOptimizationOptions', opt);

%optimization results 
optRes = Mdl.HyperparameterOptimizationResults;

best = optRes.XAtMinObjective;   % 
bestBox   = best.BoxConstraint;
bestScale = best.KernelScale;

%get RBF cv loss 
tBestRBF = templateSVM( ...
    'KernelFunction', 'rbf', ...
    'Standardize',    true, ...
    'Type',           'classification', 'BoxConstraint',bestBox, 'KernelScale',bestScale);

CVMdlRBF = fitcecoc(X_train,Y_train,'Learners',tBestRBF,  'Coding','onevsone', ...
                    'CVPartition', cv_inner);
lossRBF = kfoldLoss(CVMdlRBF);
%%lossRBF =   0.2698

%% train linear model 

tSVMlin = templateSVM( ...
    'KernelFunction', 'linear', ...
    'Standardize',    true, ...
    'Type',           'classification');
%options for automatic fitting of parameters 
opt = struct( ...
    'CVPartition', cv_inner, ...
    'Optimizer', 'bayesopt', ...
    'AcquisitionFunctionName','expected-improvement-plus', ...
    'ShowPlots', true, ...
    'Verbose', 1, ...
    'MaxObjectiveEvaluations', 30);

%parameter optimization 
[Mdl] = fitcecoc(X_train, Y_train, ...
    'Learners', tSVMlin, ...
    'Coding', 'onevsone', ...
    'OptimizeHyperparameters', {'BoxConstraint'}, ...
    'HyperparameterOptimizationOptions', opt);

%optimization results 
optRes = Mdl.HyperparameterOptimizationResults;

best = optRes.XAtMinObjective;   % 
bestBox   = best.BoxConstraint;


%get Linear model cv loss 
tBestLinear = templateSVM( ...
    'KernelFunction', 'linear', ...
    'Standardize',    true, ...
    'Type',           'classification', 'BoxConstraint',bestBox);

CVMdlLinear = fitcecoc(X_train,Y_train,'Learners',tBestLinear,  'Coding','onevsone', ...
                    'CVPartition', cv_inner);
lossLinear = kfoldLoss(CVMdlLinear);
%%lossLinear =   0.2710


%% Conclusion linear model just a bit worse than RBF model, 
%% we use LINEAT SVM model with optimized BoxConstraint 
tBest = tBestLinear;
bestParams = tBestLinear;
cvErr = lossLinear;


%% Fit and save the chosen model 


Mdl1 = fitcecoc(X_train, Y_train, 'Learners', tBest, 'Coding','onevsone', 'FitPosterior', true);

[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,400 400], 'normal', 'on');

% 1) Out-of-fold predictions
CVMdlBest = fitcecoc(X_train,Y_train,'Learners',tBest,  'Coding','onevsone', ...
                    'CVPartition', cv_inner);

yHat = kfoldPredict(CVMdlBest);   % each point predicted by a model that did not see it
%yHat = Mdl1.predict(X_train(:, :));   % each point predicted by a model that did not see it

% 2) Confusion matrix (counts)
C = confusionmat(Y_train, yHat);

% 3) Row-normalized confusion matrix (per-class mixing)
Cnorm = C ./ sum(C, 2);

% 4) Nice plot
%confusionchart(Y_train, yHat);      % by default shows counts
confusionchart(Y_train, yHat, 'Normalization','row-normalized');  % mixing per true class
title(LM_Model_title)
%Save the algorithm

save(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']), 'Mdl1', 'ClassifID', ...
    'mask_seq2', 'cv_inner', 'bestParams', 'cvErr', 'LM_Model_title');


%use Fram_0435_EnvType_PartDepend to plot partial dependencies 


%Save predicted cluster ID
tDMC.ClustID_predicted = Mdl1.predict(tStPredictors);
writenewtable(tDMC, fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID_obsrvpred.xlsx'));