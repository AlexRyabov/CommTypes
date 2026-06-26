%% Chl distribution
% Alex Ryabov, 2025

proj_Initialize

predictors = {'Lon', 'Month_doy', 'Depth'};

%Remove nans
tStPredictors = tStPredictors(~isnan(tStPredictors.chlorophyll_a_0), :);

Y = tStPredictors.chlorophyll_a_0;
X = tStPredictors(:, predictors);
%% transform the response
alpha =0.25;
var_tranf = @(y) y.^alpha ;
inv_var_tranf = @(tranf_y) tranf_y.^(1/alpha);
%[Y, lambda] = boxcox(tStPredictors.chlorophyll_a_0);
% var_tranf = @(y) log10(y);
%inv_var_tranf = @(tranf_y) 10.^tranf_y;
Y = var_tranf(tStPredictors.chlorophyll_a_0);
figure(1)
clf
histogram(Y)
%%

%train svm with automatic feature selection
% [X_train, Y_train, X_test, Y_test, cv_inner] = ...
%     ml_blockedCV(X,Y, tStPredictors.Cruise, tStPredictors.MonthNum, HGP.ML.KFolds);

cv_inner = cvpartition(tStPredictors.MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
nFolds = cv_inner.NumTestSets;

%Test and training
X_train = X(:, :);
Y_train = Y(:);
X_test = [];
Y_test = []; 

%%
LM_models = {'RF', 'SVM_auto', 'SVM', 'Boost1', 'Boost2',  'Lin_stepwise', 'GAM', 'FitNet', 'RFNative'}';
Select_Best_Model = true;

if Select_Best_Model
    %%tune the a regressor
    Models = table(LM_models);
    Models.Row = LM_models;
    %Models.LM_models = LM_models;
    for iM = 1:length(LM_models)
        iM
        [Mdl, TuneRes] = ml_tune_model(LM_models{iM},X_train,Y_train,cv_inner);
        Models.Mdl{iM} = Mdl;
        Models.TuneRes{iM} = TuneRes;
        Models.cvloss(iM) = TuneRes.bestCVLoss;
        %plot partial dependence 
        [fg, MaxFigure] = addFindFigure(['PD' LM_models{iM}], MaxFigure, [10, 20,1200 400], 'normal', 'on');
        for iP = 1:length(predictors)
            nexttile
            plotPartialDependence( Models.Mdl{iM}, predictors{iP});
        end
    end
    save(fullfile(HGP.Data.Path2ResultsData, ['trainedChlRegrModels.mat']), 'Models');
else
    tmp = load(fullfile(HGP.Data.Path2ResultsData, ['trainedChlRegrModels.mat']));
    Models = tmp.Models;
    clearvars tmp;
    BestModel = Models.Mdl{'SVM_auto'};
end
   
%% plot chl distribution 
%1 svm, manual fit
%4 boost_2, ext to  EGC
%5 RF, max in the middle

iM = 3;
BestModel = Models.Mdl{iM};
Models.LM_models{iM}


[fg, MaxFigure] = addFindFigure('Chl_distr', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');

depths = [0, 30; 30, 60; 60, 100];
titles = {'Shallow (0–30 m)', 'Middle (30-60m)', 'Deep (60–100 m)'};
Lon0 = 0;
Month0 = 7;
lonVec   = linspace(min(tStPredictors.Lon), max(tStPredictors.Lon), 50)';     %
monthVec   = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 50)';     %
[LonGrid, MonthGrid] = meshgrid(lonVec, monthVec);

fun = @(Tbl) predict(BestModel, Tbl);

pds = cell(height(depths));
for iD = 1:height(depths)

    depthGrid = (depths(iD, 1):5:depths(iD, 2))';
    selectedData = table( ...
        repmat(Lon0,   size(depthGrid)), ...
        repmat(Month0, size(depthGrid)), ...
        depthGrid, ...
        'VariableNames', {'Lon','Month_doy','Depth'});

    [pd,x,y] = partialDependence(fun, {'Lon', 'Month_doy'}, selectedData, ...
        "QueryPoints", [LonGrid(:) MonthGrid(:)], ...
        "UseParallel",true);
    pd = inv_var_tranf(pd);
    pds{iD} = pd;

    nexttile(iD)
    
    contourf(x, y, pd, 10, 'ShowText','on');      % x -> Lon, y -> Month_doy

    %if iD == height(depths)
        colorbar;
    %end

    xlabel('Longitude');
    ylabel('Month');
    title(titles(iD));
    % hold on;
    % contour(x, y, pMax2D, [0, 0.4, 0.6, 0.8, 1], 'k', 'LineWidth', 1, 'ShowText','on');
    % hold off;
end

   %plot partial dependence 
        [fg, MaxFigure] = addFindFigure(['PD' Models.LM_models{iM}], MaxFigure, [10, 20,1200 400], 'normal', 'on');
        for iP = 1:length(predictors)
            nexttile
            plotPartialDependence( Models.Mdl{iM}, predictors{iP});
        end


%% plot chl from data 
[fg, MaxFigure] = addFindFigure('Chl_data', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');

depths = [0, 30; 30, 60; 60, 100];
titles = {'Shallow (0–30 m)', 'Middle (30-60m)', 'Deep (60–100 m)'};

pds = cell(height(depths));
for iD = 1:height(depths)
    ind = tStPredictors.Depth >= depths(iD, 1) & tStPredictors.Depth <= depths(iD, 2);
    selectedData =   tStPredictors(ind, {predictors{:}, 'chlorophyll_a_0'});
  
   nexttile
   scatter(selectedData.Lon, selectedData.Month_doy, 30, selectedData.chlorophyll_a_0,  'filled');   

   caxis([0, 2])
        colorbar;

    xlabel('Longitude');
    ylabel('Month');
    title(titles(iD));
    % hold on;
    % contour(x, y, pMax2D, [0, 0.4, 0.6, 0.8, 1], 'k', 'LineWidth', 1, 'ShowText','on');
    % hold off;
end



