%Compare chl prediction with and without community type
%Use Month_doy, Lon, temp_imputed, depth  + comm type as a parameter 

%Train one model which includes only environment and another with Type ID
%as a random parameters 
%do it for Month_doy, Lon, temp_imputed, depth  and Month_doy, temp_imputed, depth 

% Alex Ryabov, 2025

proj_Initialize

predictors = {'Lon', 'Month_doy', 'Depth'};

%Load diffusion map clusters
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID_obsrvpred.xlsx'));
tStPredictors.ClustID = categorical(tDMC.ClustID);
%%mean Depth =  41.5035
tStPredictors.Depth = tStPredictors.Depth - mean(tStPredictors.Depth);
tStPredictors.ClustID_predicted = categorical(tDMC.ClustID_predicted);

%Remove nans
tStPredictors = tStPredictors(~isnan(tStPredictors.chlorophyll_a_0), :);


%% transform the response
alpha =0.25;
%alpha =1;
var_tranf = @(y) y.^alpha ;
inv_var_tranf = @(tranf_y) tranf_y.^(1/alpha);
%[Y, lambda] = boxcox(tStPredictors.chlorophyll_a_0);
% var_tranf = @(y) log10(y);
%inv_var_tranf = @(tranf_y) 10.^tranf_y;
tStPredictors.transfChl = var_tranf(tStPredictors.chlorophyll_a_0);
figure(1)
clf
histogram(tStPredictors.transfChl)


%% mixed model for chl
%make and clean input table
%% --- Fit mixed model (run once) ---

% Fixed-effect-only model for marginal line
lme0 = fitlme(tStPredictors, 'transfChl ~ Depth+  MonthNum');  %%Dependence on Depth^2 is insignificant
% Mixed model: random intercept only based on predicted class
stdNoType = std(lme0.predict(tStPredictors)-tStPredictors.transfChl)

lme1Pred = fitlme(tStPredictors, 'transfChl ~ Depth +  MonthNum + (1|ClustID_predicted)');
lme2Pred = fitlme(tStPredictors, 'transfChl ~ Depth +  MonthNum + (Depth|ClustID_predicted)');
compare(lme0, lme1Pred)
compare(lme1Pred, lme2Pred)
stdWType = std(lme1Pred.predict(tStPredictors)-tStPredictors.transfChl)
%% Temperature
lme0_temp = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+temp_imputed', 'FitMethod','REML')  %%Dependence on Depth^2 is insignificant
% Mixed model: random intercept only based on predicted class
stdNoType_temp = std(lme0_temp.predict(tStPredictors)-tStPredictors.transfChl)

%lme1Pred_temp = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+temp_imputed+(Depth|ClustID_predicted)', 'FitMethod','REML')

%compare(lme0_temp, lme1Pred_temp)
%stdWType_temp = std(lme1Pred_temp.predict(tStPredictors)-tStPredictors.transfChl)

lme1Pred_temp_nocor = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+temp_imputed+(1|ClustID_predicted)+(Depth-1|ClustID_predicted)', 'FitMethod','REML')
stdWType_temp_nocorr = std(lme1Pred_temp_nocor.predict(tStPredictors)-tStPredictors.transfChl)
compare(lme0_temp, lme1Pred_temp_nocor)

%true community type 
lme1Pred_temp_nocor_true = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+temp_imputed+(1|ClustID)+(Depth-1|ClustID)', 'FitMethod','REML')
stdWType_temp_nocorr_true = std(lme1Pred_temp_nocor_true.predict(tStPredictors)-tStPredictors.transfChl)
compare(lme0_temp, lme1Pred_temp_nocor_true)


lme1Pred_temp_nocor_true_norand = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+temp_imputed+ClustID+ Depth:ClustID', 'FitMethod','ML')
stdWType_temp_nocorr_true = std(lme1Pred_temp_nocor_true.predict(tStPredictors)-tStPredictors.transfChl)
compare(lme0_temp, lme1Pred_temp_nocor_true)


%RE = randomEffects(lme1Pred_temp)
%% test with leave one year out 
years = unique(tStPredictors.Year);
rmse0 = nan(numel(years),1);
rmse1 = nan(numel(years),1);

form0 = 'transfChl ~ Depth + MonthNum + temp_imputed';
form1 = 'transfChl ~ Depth + MonthNum + temp_imputed + (1|ClustID) + (Depth-1|ClustID)';

for i = 1:numel(years)
    y = years(i);
    tr = tStPredictors.Year ~= y;
    te = tStPredictors.Year == y;

    % Fit on training only
    m0 = fitlme(tStPredictors(tr,:), form0, 'FitMethod','REML');
    m1 = fitlme(tStPredictors(tr,:), form1, 'FitMethod','REML');

    % Predict on held-out year
    yte = (tStPredictors.transfChl(te));
    yhat0 = (predict(m0, tStPredictors(te,:)));
    yhat1 = (predict(m1, tStPredictors(te,:)));

    rmse0(i) = sqrt(mean((yhat0 - yte).^2));
    rmse1(i) = sqrt(mean((yhat1 - yte).^2));
end

meanRMSE0 = mean(rmse0)
meanRMSE1 = mean(rmse1)
deltaRMSE = meanRMSE0 - meanRMSE1   % >0 means CommType helps

delta = rmse0 - rmse1;  % improvement per held-out year

meanDelta = mean(delta);
medDelta  = median(delta);
nPos      = sum(delta > 0);
nTot      = numel(delta);

fprintf("Mean ΔRMSE = %.4f\n", meanDelta);
fprintf("Median ΔRMSE = %.4f\n", medDelta);
fprintf("Years improved: %d / %d\n", nPos, nTot);

B = 5000;
bootMean = nan(B,1);
n = numel(delta);

for b = 1:B
    idx = randi(n, n, 1);      % resample years with replacement
    bootMean(b) = mean(delta(idx));
end

CI = prctile(bootMean, [2.5 97.5]);
fprintf("Bootstrap 95%% CI for mean ΔRMSE: [%.4f, %.4f]\n", CI(1), CI(2));

p = signrank(delta, 0, 'tail', 'right');  % H1: median(delta) > 0
fprintf("Wilcoxon signed-rank p-value (ΔRMSE>0): %.4g\n", p);

rel = deltaRMSE / meanRMSE0;
fprintf("Relative RMSE improvement: %.2f%%\n", 100*rel);

%% Longitude and Lat do not give a significant effect 
lme0_temp_lon = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+Lon+temp_imputed', 'FitMethod','REML')  %%
lme0_temp_lat = fitlme(tStPredictors, 'transfChl ~ Depth+MonthNum+Lat+temp_imputed', 'FitMethod','REML')  %%
compare(lme0_temp, lme1Pred_temp_nocor)

%% --- Plot --- plot works only for depth
if 1
    [fg, MaxFigure] = addFindFigure('Chl_rand_effects', MaxFigure, [10, 20, 1.1003e+03 310.6667], 'normal', 'on');

%% Random effects visualization using Stat output from randomEffects (MATLAB R2025b)
% Required variables:
%   - lme1Pred_temp_nocor : fitted LinearMixedModel
%   - tStPredictors       : table used for fitting (same variables as in the model)
%
% This script visualizes:
%   (1) Random intercept shifts (BLUPs) with 95% CI
%   (2) Random Depth slope shifts (BLUPs) with 95% CI
%   (3) Cluster-specific intercept and Depth slope = fixed + random
%   (4) Cluster-specific prediction lines vs Depth (MonthNum and temp fixed)

%lme = lme1Pred_temp_nocor_true;
lme = lme1Pred_temp_nocor_true;
tbl = tStPredictors;

slopeVar = 'Depth';

%%1) Get random effects stats (Stat contains Level, Name, Estimate, Lower, Upper, etc.)
[RandEff, Names, Stat] = randomEffects(lme, 'Alpha', 0.05); 

% Stat is a dataset array (or table-like object) with fields:
% Group, Level, Name, Estimate, SEPred, tStat, DF, pValue, Lower, Upper

%%2) Extract rows for (Intercept) and for Depth
isInt   = strcmp(Stat.Name, '(Intercept)');
isSlope = strcmp(Stat.Name, slopeVar);

S0 = Stat(isInt, :);
SD = Stat(isSlope, :);

% Convert cluster levels to strings for labeling
cl0 = string(S0.Level);
clD = string(SD.Level);

% Align slope rows to intercept rows by cluster level
[tf, idx] = ismember(cl0, clD);
if any(~tf)
    error('Some cluster levels are missing in the slope random effects.');
end
SD = SD(idx, :);

% Build a compact structure/table for plotting
Clust = cl0;
b0    = S0.Estimate;
L0    = S0.Lower;
U0    = S0.Upper;

bD    = SD.Estimate;
LD    = SD.Lower;
UD    = SD.Upper;

p0    = S0.pValue;
pD    = SD.pValue;

%%3) Caterpillar plots with 95% CI (from Lower/Upper)
clf
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% --- Random intercept shifts ---
nexttile(1);
[~,ord0] = sort(b0,'descend');
ord0 = 1:numel(bD);
y = 1:numel(b0);

x  = b0(ord0);
xl = x - L0(ord0);   % left error (Estimate - Lower)
xu = U0(ord0) - x;   % right error (Upper - Estimate)

errorbar(x, y, xl, xu, 'horizontal','o', ...
    'MarkerSize',6,'MarkerFaceColor',[0.2 0.45 0.9], ...
    'Color',[0.2 0.45 0.9],'LineWidth',1.2);
xline(0,'k--','LineWidth',1);
grid on;
yticks(y);
yticklabels(Clust(ord0));
xlabel('Random intercept shift (BLUP)');
title('Random intercept shifts by cluster (95% CI)');
ylim([0.5, numel(b0)+0.5])

% --- Random slope shifts (Depth) ---
nexttile(2);
[~,ordD] = sort(bD,'descend');
ordD = 1:numel(bD);
y = 1:numel(bD);

x  = bD(ordD);
xl = x - LD(ordD);
xu = UD(ordD) - x;

errorbar(x, y, xl, xu, 'horizontal','o', ...
    'MarkerSize',6,'MarkerFaceColor',[0.15 0.65 0.35], ...
    'Color',[0.15 0.65 0.35],'LineWidth',1.2);
xline(0,'k--','LineWidth',1);
grid on;
yticks(y);
yticklabels(Clust(ordD));
xlabel('Random slope shift for Depth (BLUP)');
title('Random Depth slope shifts by cluster (95% CI)');
ylim([0.5, numel(bD)+0.5])

%% 4) Cluster-specific intercept and slope = fixed + random
beta = fixedEffects(lme);
cn   = string(lme.CoefficientNames);

b0fix = beta(cn == "(Intercept)");
bMfix = beta(cn == "MonthNum");
bDfix = beta(cn == slopeVar);
bTfix = beta(cn == "temp_imputed");


Intercept_cluster = b0fix + b0;
SlopeDepth_cluster = bDfix + bD;



nexttile(3);
[~,o1] = sort(Intercept_cluster,'descend');
barh(categorical(Clust(o1)), Intercept_cluster(o1), 'FaceColor',[0.2 0.45 0.9]);
grid on;
xlabel('Intercept (fixed + random)');
title('Cluster-specific intercept');

nexttile(4);
[~,o2] = sort(SlopeDepth_cluster,'descend');
barh(categorical(Clust(o2)), SlopeDepth_cluster(o2), 'FaceColor',[0.15 0.65 0.35]);
grid on;
xlabel('Depth slope (fixed + random)');
title('Cluster-specific Depth slope');

%% 5) Cluster-specific prediction lines vs Depth (MonthNum and temp fixed)
mn  = mean(tbl.MonthNum, 'omitnan');
tmp = mean(tbl.temp_imputed, 'omitnan');

dmin = min(tbl.Depth);
dmax = max(tbl.Depth);
dgrid = linspace(dmin, dmax, 80)';

figure('Color','w'); hold on; grid on;

% Fixed-effects line
y_fixed = b0fix + bMfix*mn + bTfix*tmp + bDfix*dgrid;
plot(dgrid, y_fixed, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Fixed effects');

% Cluster-specific lines
cols = lines(numel(Clust));
for i = 1:numel(Clust)
    y_i = (Intercept_cluster(i) + bMfix*mn + bTfix*tmp) + SlopeDepth_cluster(i)*dgrid;
    plot(dgrid, y_i, '-', 'Color', cols(i,:), 'LineWidth', 1.6, ...
        'DisplayName', "Cluster " + Clust(i));
end

xlabel('Depth');
ylabel('Predicted transfChl');
title('Depth \rightarrow transfChl: cluster-specific lines (MonthNum,temp fixed)');
legend('Location','bestoutside');

%% 6) Optional: annotate clusters with p-values (random effects "t-tests")
% Note: p-values for BLUPs are often not the main inferential target in mixed models.
% Uncomment if you still want to display them in the command window.
%{
disp(table(Clust, b0, p0, bD, pD, ...
    'VariableNames', {'Cluster','RE_Intercept','p_Intercept','RE_Depth','p_Depth'}));
%}

%% 7) Optional: refit with centered Depth for more interpretable intercepts
% Centering makes the intercept correspond to "average depth" rather than Depth=0.
%{
tbl.Depth_c = (tbl.Depth - mean(tbl.Depth,'omitnan'));
lme_c = fitlme(tbl, ...
    'transfChl ~ Depth_c + MonthNum + temp_imputed + (1|ClustID_predicted) + (Depth_c-1|ClustID_predicted)', ...
    'FitMethod','REML');
% Then rerun this script with slopeVar = 'Depth_c' and lme = lme_c.
%}
    
end
%%

