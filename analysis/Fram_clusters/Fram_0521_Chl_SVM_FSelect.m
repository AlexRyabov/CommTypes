%% Chl distribution as a function of Lon, Month, and Depth and other env variables
%with feature selection of the best combination
%feature selection is not so stable and you need to repeat is few times to
%get the final results
% Alex Ryabov, 2025
clearvars
proj_Initialize

%Load diffusion map clusters
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
tStPredictors.ClustID = categorical(tDMC.ClustID);

% Model ID

modelID = 'SimpleRobust'; %for feature selection, blocked CV
%modelID = 'SmoothTuning'; %Use selected features, random CV

switch modelID
    case 'SimpleRobust'
        predictors = {'Month_doy', 'Depth', 'temp_imputed', 'ClustID', 'sea_ice_now', 'sea_ice_monthly', 'nitrate_imputed'};
        %train svm with automatic feature selection
        XY =  rmmissing(tStPredictors(:, {predictors{:}, 'Cruise', 'MonthNum',  'chlorophyll_a_0'}));
        [foldIdx, cruiseFold] = makeCruiseMonthFolds(XY.Cruise, XY.MonthNum, HGP.ML.KFolds);
        cv_inner = cvpartition("CustomPartition", foldIdx);

    case 'SmoothTuning'
        predictors = {'Month_doy', 'Depth', 'temp_imputed', 'ClustID'};
        XY =  rmmissing(tStPredictors(:, {predictors{:}, 'Cruise', 'MonthNum',  'chlorophyll_a_0'}));
        cv_inner = cvpartition(XY.MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
    otherwise
        error('')
end

Y = XY.chlorophyll_a_0;
X = XY(:, predictors);

%Remove nans
%%transform the response
alpha =0.2;
%alpha =1;
var_tranf = @(y) y.^alpha ;
inv_var_tranf = @(tranf_y) tranf_y.^(1/alpha);

Y = var_tranf(XY.chlorophyll_a_0);
[fg, MaxFigure] = addFindFigure('Chl_rand_effects', MaxFigure, [10, 20, 400 400]);

clf
histogram(Y)







%% tune the model
tuneModel = true;
LM_model = 'SVM';
rng(1);
%%tune model parameter
[BestModel, TuneRes] = ml_tune_model(LM_model,X_train,Y_train,cv_inner);

%%
% bp.BoxConstraint = 100;
% bp.KernelScale = 50;
bp = TuneRes.bestPoint;

%training function
ml_train_func = @(X_tr, Y_tr, bestPoint) fitrsvm( ...
    X_tr, Y_tr, ...
    'Standardize', true, ...                  % standardize predictors
    'KernelFunction', 'gaussian', ...         % or 'linear', 'rbf' etc.
    'BoxConstraint',  bestPoint.BoxConstraint,...
    'KernelScale',  bestPoint.KernelScale);

%criterion for feature selection (cv losss)
fsCritFun = @(X_tr, Y_tr,X_te, Y_te) ...
    mlTrainAndLoss(X_tr, Y_tr, X_te, Y_te, ml_train_func, bp);
if 0  %use dummy variable for Clust_ID and use it for feature selection
    %When clust ID is involved into selection, the model is very
    %complicated and feature selection does not work
    X_trainClustID  =  onehotencode(X_train(:, 'ClustID'));
    X_trainNew = X_train;
    X_trainNew.ClustID = [];
    X_trainNew = [X_trainNew, X_trainClustID];

    %X_train.ClustID = double(X_train.ClustID);
    [mask_seq, history] = sequentialfs(fsCritFun, X_trainNew{:, :}, Y_train, 'CV', cv_inner,"Direction", "backward", "Options",statset("Display","iter", "UseParallel",true));

    clNames = arrayfun(@num2str, 1:HGP.dmCluster.NumComponents, 'UniformOutput', false);
    X_trainNew.Properties.VariableNames(mask_seq)
    predictors_sel =   X_trainNew.Properties.VariableNames(mask_seq);
    [predictors_sel, ia] = setdiff(predictors_sel, clNames, 'stable');
    if length(predictors_sel) < length(X_trainNew.Properties.VariableNames(mask_seq))
        predictors_sel =  {predictors_sel{:}, 'ClustID'};
    end
    [BestModel, TuneRes] = ml_tune_model(LM_model,X_train(:, predictors_sel),Y_train,cv_inner);


    %  {'Month_doy'}    {'Depth'}    {'temp_imputed'}    {'sea_ice_now'}    {'sea_ice_monthly'}    {'ClustID'}
    %   bestCVLoss: 0.0251


    X_trainClustID  =  onehotencode(X_train(:, 'ClustID'));
    X_trainNew = X_train(:, setdiff(predictors_sel, 'ClustID', 'stable'));
    X_trainNew = [X_trainNew, X_trainClustID];

    [mask_seq, history] = sequentialfs(fsCritFun, X_trainNew{:, :}, Y_train, 'CV', cv_inner,"Direction", "backward", "Options",statset("Display","iter", "UseParallel",true));

    clNames = arrayfun(@num2str, 1:HGP.dmCluster.NumComponents, 'UniformOutput', false);
    X_trainNew.Properties.VariableNames(mask_seq)
    predictors_sel =   X_trainNew.Properties.VariableNames(mask_seq);
    [predictors_sel, ia] = setdiff(predictors_sel, clNames, 'stable');
    if length(predictors_sel) < length(X_trainNew.Properties.VariableNames(mask_seq))
        predictors_sel =  {predictors_sel{:}, 'ClustID'};
    end
    [BestModel, TuneRes] = ml_tune_model(LM_model,X_train(:, predictors_sel),Y_train,cv_inner);
    %bestCVLoss: 0.0250
    %  {'Month_doy'}    {'Depth'}    {'temp_imputed'}    {'sea_ice_now'}    {'sea_ice_monthly'}    {'ClustID'}
else
    %We do not use cluster ID, just select from environmental predictors
    %finally we end up with D, T, Doy
    canRemove = true;
    [predictors_sel, ia] = setdiff(predictors, 'ClustID', 'stable');
    while canRemove
        [BestModel, TuneRes] = ml_tune_model(LM_model,X_train(:, predictors_sel),Y_train,cv_inner);
        fsCritFun = @(X_tr, Y_tr,X_te, Y_te) ...
            mlTrainAndLoss(X_tr, Y_tr, X_te, Y_te, ml_train_func, TuneRes.bestPoint);
        X_trainNew  = X_train(:, predictors_sel);
        [mask_seq, history] = sequentialfs(fsCritFun, X_trainNew{:, :}, Y_train, 'CV', cv_inner,"Direction", "backward", "Options",statset("Display","iter", "UseParallel",true));
        canRemove = sum(mask_seq) < numel(predictors_sel);
        predictors_sel = X_trainNew.Properties.VariableNames(mask_seq)
    end
    %%
    [BestModel, TuneRes] = ml_tune_model(LM_model,X_train(:, predictors_sel),Y_train,cv_inner);
    % blocked CV
    % bestCVLoss: 0.0252,   {'Month_doy'}    {'Depth'}    {'temp_imputed'}
    % random CV
    %  bestCVLoss: 0.0198   {'Month_doy'}    {'Depth'}    {'temp_imputed'}

    %Selected features + cluster ID
    [BestModelClust, TuneResClust] = ml_tune_model(LM_model,X_train(:, [predictors_sel, "ClustID"]),Y_train,cv_inner);
    % blocked CV
    %            bestCVLoss: 0.0248
    % random CV
    % bestCVLoss: 0.0183
end



%%
[fg, MaxFigure] = addFindFigure(['predVSObs'], MaxFigure, [10, 20,1200 400]);
%plot predicted versus observed values of chl
%yHat = BestModel.predict(X_train);
mdlCV = crossval(BestModel, 'CVPartition', cv_inner);
yHat  = kfoldPredict(mdlCV);

%inv_trans2plot = @(X) inv_var_tranf(X);
inv_trans2plot = @(X) X;

chlObs = inv_trans2plot(Y_train);
yHat = inv_trans2plot(yHat);


nexttile
scatter(yHat, chlObs)
hold on
hold on;
plot([min(yHat), max(yHat)], [min(yHat), max(yHat)], 'r--'); % 45-degree line
hold off;
ylabel('Observed chl-a')
xlabel('Predicted chl-a')

t2 = nexttile(2)
residuals = chlObs -yHat ;
scatter(yHat, residuals, 30,  -X_train.Depth, 'filled');
xlabel('Predicted chl-a');
ylabel('Residual (Obs -Pred)');
%colormap(t2, "jet")
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
cb = colorbar
cb.Label.String = 'Depth (m)';


t2 = nexttile(3)
clrs = f_Clrs_fresh(HGP.dmCluster.NumComponents);

vs = boxplot(residuals, X_train.ClustID);
%yscale('log')
ylabel('Residual (Obs -Pred)')
xlabel('Community type')

%% plot partial dependence
[fg, MaxFigure] = addFindFigure(['PD' Models.LM_models{iM}], MaxFigure, [10, 20,1200 400], 'normal', 'on');
for iP = 1:length(predictors_sel)
    nexttile
    plotPartialDependence(BestModel, predictors_sel{iP});
end







