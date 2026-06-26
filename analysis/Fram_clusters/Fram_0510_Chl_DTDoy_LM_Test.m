%% Chl distribution as a function of Lon, Month, and Depth, SVM fit and residiaul by type
% Alex Ryabov, 2025

proj_Initialize

predictors = {'Lon', 'Month_doy', 'Depth'};

%Load diffusion map clusters
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
tStPredictors.ClustID = categorical(tDMC.ClustID);



%Remove nans
tStPredictors = tStPredictors(~isnan(tStPredictors.chlorophyll_a_0), :);


Y = tStPredictors.chlorophyll_a_0;
X = tStPredictors(:, predictors);
%% transform the response
alpha =0.25;
%alpha =1;
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
LM_models = {'RF', 'SVM_auto', 'SVM'}'; %, 'Boost1', 'Boost2',  'Lin_stepwise', 'GAM', 'FitNet', 'RFNative'}';
Select_SaveBest_Model = false;

if Select_SaveBest_Model
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

iM = 3;  %
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
    if iD == 1  %plot ice edge line

        x = tIceEdgeLoc.ICE20LongMin;
        y = tIceEdgeLoc.Month;
        x = [x; tIceEdgeLoc.ICE20LongMax(end:-1:1)];
        y = [y; tIceEdgeLoc.Month(end:-1:1)];
        p = patch(x', y'+0.5, 'r',  'EdgeColor', 'b');
        p.LineWidth = 3;
        p.FaceAlpha = 0.0;
    end     
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
    selectedData =   tStPredictors(ind, [predictors(:)', {'chlorophyll_a_0'}]);

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


%% plot predicted vs observed data
[fg, MaxFigure] = addFindFigure('Chl_data_res', MaxFigure, [-1, -1, 1261 377], 'normal', 'on');


chlObs = (tStPredictors.chlorophyll_a_0);
chlHat = inv_var_tranf(fun(tStPredictors));


% Plot predicted vs observed data
t = nexttile
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
residuals = var_tranf(tStPredictors.chlorophyll_a_0) -(chlHat)  ;
scatter(tStPredictors.chlorophyll_a_0, residuals, 30,  -tStPredictors.Depth, 'filled');
xlabel('Observed chlorophyll-a (transformed)');
ylabel('Residual (Obs -Pred)');
%colormap(t2, "jet")
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
cb = colorbar
cb.Label.String = 'Depth (m)';
%%
%
%

% Calculate the residuals


residuals = (tStPredictors.chlorophyll_a_0) - inv_var_tranf(chlHat);
%residuals = var_tranf(tStPredictors.chlorophyll_a_0) - (chlHat);
% Pairwise comparisons (Tukey-Kramer)

nexttile(3)
[p, tblA, stats] = anova1(residuals, tStPredictors.ClustID, 'off');  % 'on' shows built-in plot
c = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
% c = [g1 g2 lowerCI meanDiff upperCI pValue]
disp(p);
cMulty = multcompare(stats, 'CType','tukey-kramer', "Alpha" ,0.05 );
sig = cMulty(cMulty(:,6) < 0.05, :)


p_levene = vartestn(residuals, tStPredictors.ClustID, 'TestType','LeveneAbsolute', 'Display','off');
p_levene



p_kw = kruskalwallis(residuals, categorical(tStPredictors.ClustID), 'off');
p_kw

Type = categorical(tStPredictors.ClustID);
tbl = table(residuals, Type, 'VariableNames', {'Resid','Type'});
lm = fitlm(tbl, 'Resid ~ Type', 'RobustOpts','on');
anova(lm,'summary')



%%
nexttile(3)
[p, tblA, stats] = anova1(residuals, tStPredictors.ClustID, 'off');  % 'on' shows built-in plot
c = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
% c = [g1 g2 lowerCI meanDiff upperCI pValue]
disp(p);
cMulty = multcompare(stats, 'CType','tukey-kramer', "Alpha" ,0.05 );
sig = cMulty(cMulty(:,6) < 0.05, :)


p_levene = vartestn(residuals, tStPredictors.ClustID, 'TestType','LeveneAbsolute', 'Display','off');
p_levene





p_kw = kruskalwallis(residuals, categorical(tStPredictors.ClustID), 'off');
p_kw

Type = categorical(tStPredictors.ClustID);
tbl = table(residuals, Type, 'VariableNames', {'Resid','Type'});
lm = fitlm(tbl, 'Resid ~ Type', 'RobustOpts','on');
anova(lm,'summary')


%%%
%% residuals vs water type
[fg, MaxFigure] = addFindFigure('Chl_res_wtype', MaxFigure, [-1, -1, 651 564], 'normal', 'on');
t = nexttile
tStPredictors.wTypec = categorical(tStPredictors.wType,  {'PSW', 'wPSW', 'MAW', 'AW'}, {'Polar w', 'warm Pol w', 'mod Atl w', 'Atl w'});



j = (rand(height(tStPredictors), 2)-0.5)*0.7;
s = scatter(double(tStPredictors.ClustID) + j(:, 1), double(tStPredictors.wTypec) + j(:, 2), 30, residuals, 'filled')
s.MarkerFaceAlpha = 0.7;
colormap(redblue)
colorbar
caxis([-2, 2])
ax = gca
ax.XTick = [1:8]
ax.YTick = [1:4]
ax.YTickLabel = categories(tStPredictors.wTypec);
ax.YTickLabelRotation(:) = 30;
xlabel('Community type')
ylabel('Water type')
view([90, 90])


