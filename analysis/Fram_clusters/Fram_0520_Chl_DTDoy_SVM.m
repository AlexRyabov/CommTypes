%% Chl distribution as a function of Lon, Month, and Depth
% Alex Ryabov, 2025

proj_Initialize

predictors = {'Lon', 'Month_doy', 'Depth'};

%Load diffusion map clusters
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
tStPredictors.ClustID = categorical(tDMC.ClustID);



%Remove nans
tStPredictors = tStPredictors(~isnan(tStPredictors.chlorophyll_a_0), :);


Y = tStPredictors.chlorophyll_a_0;
X = tStPredictors(:, {predictors{:}, 'ClustID'});
%% transform the response
alpha =0.25;
%alpha =1;
var_tranf = @(y) y.^alpha ;
inv_var_tranf = @(tranf_y) tranf_y.^(1/alpha);
%[Y, lambda] = boxcox(tStPredictors.chlorophyll_a_0);
% var_tranf = @(y) log10(y);
%inv_var_tranf = @(tranf_y) 10.^tranf_y;
Y = var_tranf(tStPredictors.chlorophyll_a_0);
 [fg, MaxFigure] = addFindFigure('Chl_rand_effects', MaxFigure, [10, 20, 400 400]);

clf
histogram(Y)
%%
%set MonthNum to 6 for end of May (4 observations)
tStPredictors.MonthNum(tStPredictors.Month_doy<6) = round(tStPredictors.Month_doy(tStPredictors.Month_doy<6));

%train svm with automatic feature selection
 [X_train, Y_train, X_test, Y_test, cv_inner] = ...
ml_blockedCV(X,Y, tStPredictors.Cruise, tStPredictors.MonthNum, HGP.ML.KFolds);


%cv_inner = cvpartition(tStPredictors.MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
nFolds = cv_inner.NumTestSets;

%Test and training
X_train = X(:, :);
Y_train = Y(:);
X_test = [];
Y_test = [];




%%
tuneModel = true;
%LM_model = 'SVMrbf';
LM_model = 'SVMlinear';


if tuneModel
    %%tune the a regressor
    rng(1);
    [BestModel, TuneRes] = ml_tune_model(LM_model,X_train(:, predictors),Y_train,cv_inner);
    [fg, MaxFigure] = addFindFigure(['PD' LM_model], MaxFigure, [10, 20,1200 400], 'normal', 'on');
    for iP = 1:length(predictors)
        nexttile
        plotPartialDependence( BestModel, predictors{iP});
    end
    save(fullfile(HGP.Data.Path2ResultsData, ['trainedChlDTDoyBest.mat']), 'BestModel', 'TuneRes', 'LM_model');
else
    tmp = load(fullfile(HGP.Data.Path2ResultsData, ['trainedChlDTDoyBest.mat']));
    BestModel = tmp.BestModel;
    TuneRes = tmp.TuneRes;
    LM_model = tmp.LM_model;
    clearvars tmp;
end

%% plot chl distribution by depth, doy and Lon

[fg, MaxFigure] = addFindFigure('Chl_distr', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');

depths = [0, 30; 30, 60; 60, 100];
titles = {'Shallow (0–30 m)', 'Middle (30-60m)', 'Deep (60–100 m)'};
Lon0 = 0;
Month0 = 7;
lonVec   = linspace(min(tStPredictors.Lon), max(tStPredictors.Lon), 50)';     %
monthVec   = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 50)';     %
[LonGrid, MonthGrid] = meshgrid(lonVec, monthVec);

pred_fun = @(Tbl) predict(BestModel, Tbl);

pds = cell(height(depths));
for iD = 1:height(depths)

    depthGrid = (depths(iD, 1):5:depths(iD, 2))';
    selectedData = table( ...
        repmat(Lon0,   size(depthGrid)), ...
        repmat(Month0, size(depthGrid)), ...
        depthGrid, ...
        'VariableNames', {'Lon','Month_doy','Depth'});

    [pd,x,y] = partialDependence(pred_fun, {'Lon', 'Month_doy'}, selectedData, ...
        "QueryPoints", [LonGrid(:) MonthGrid(:)], ...
        "UseParallel",true);
    pd = inv_var_tranf(pd);
    pds{iD} = pd;

    nexttile(iD)

    contourf(x, y, pd, [0:0.25:2.25], 'ShowText','on');      % x -> Lon, y -> Month_doy

    %if iD == height(depths)
    colorbar;
    %end
    caxis([0, 2.2])
    xlabel('Longitude');
    ylabel('Month');
    title(titles(iD));
    % hold on;
    % contour(x, y, pMax2D, [0, 0.4, 0.6, 0.8, 1], 'k', 'LineWidth', 1, 'ShowText','on');
    % hold off;
    
end


%% plot predicted vs observed data
mdlCV = crossval(BestModel, 'CVPartition', cv_inner);
yhat  = kfoldPredict(mdlCV);

[fg, MaxFigure] = addFindFigure('Chl_data_res', MaxFigure, [-1, -1, 1261 377], 'normal', 'on');

%inv_trans2plot = @(X) inv_var_tranf(X);
inv_trans2plot = @(X) X;

chlObs = inv_trans2plot(Y_train);
chlHat = inv_trans2plot(yhat);

% Plot predicted vs observed data
t = nexttile(1)
%colormap(t, f_Clrs_fresh(HGP.dmCluster.NumComponents));
scatter(chlObs, chlHat, 30,  -tStPredictors.Depth, 'filled');
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
hold on;
plot([min(chlObs), max(chlObs)], [min(chlObs), max(chlObs)], 'r--'); % 45-degree line
hold off;
xlabel('Observed chlorophyll-a ');
ylabel('Predicted chlorophyll-a ');
f_Set_scales({'logx', 'logy'})
%axis equal;
%grid on;
axis([0, 5, 0, 5])
cb = colorbar
cb.Label.String = 'Depth (m)';

t2 = nexttile(2)
residuals = chlObs -chlHat ;
rel_residuals = (chlObs -chlHat)./((chlObs + chlHat )*0.5);
scatter(chlObs, residuals, 30,  -tStPredictors.Depth, 'filled');
xlabel('Observed chlorophyll-a');
ylabel('Residual (Obs -Pred)');
%colormap(t2, "jet")
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
cb = colorbar
cb.Label.String = 'Depth (m)';

%%


%%
[fg, MaxFigure] = addFindFigure('Chl_data_res', MaxFigure, [-1, -1, 1261 377], 'normal', 'on');
clrs = f_Clrs_fresh(HGP.dmCluster.NumComponents);

vs = violinplot(residuals, X_train.ClustID, 'ViolinColor', clrs);
%yscale('log')
ylabel('Chl-a')
xlabel('Community type')
%%%

