%predict chl based on the probability of community type
proj_Initialize;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
%LM_Model = 'Boost_all'; %boosted regression with all predictors
%LM_Model = 'Boost_noimp'; %boosted regression no imputed predict-
%LM_Model = 'Boost_loc'; %use date and location
%LM_Model = 'Boost_env'; %use date and location
%LM_Model = 'Boost_env_sat'; %CTD, nutrients, satellite
%LM_Model = 'Boost_regular'; %regularly sampled data and most important (based on boosted )
%LM_Model = 'Boost_regular_no_depthmonth'; %regularly sampled data and most important, without depth and month info
%LM_Model = 'Boost_ctd_sat'; %regularly sampled ctd and satellite data
%LM_Model = 'Boost_ctd_sat_no_depth,doy' %regularly sampled ctd and satellite data, no depth and doy
%LM_Model = 'Boost_sat_depth'; %satellite data and depth
%LM_Model = 'Boost_selected'; %start from regular sampled and then remove less important
%LM_Model = 'Boost_most'; %satellite, imputed nutrients, depth
%LM_Model = 'Boost_best_varselect';  %after variable selection
%LM_Model = 'Boost_best_varselect+depth'  %after variable selection + depth
%LM_Model = 'Temp,Nutr,depth,Doy'  %Temp, nutr, depth and doy
%LM_Model = 'Boost_RegAndImputed'; LM_Model_title = 'Env, depth, doy'; %regularly sampled ice data, imputed environment, depth and doy
LM_Model = 'Boost_RegAndImputedNew'; LM_Model_title = 'Env, depth, doy'; %regularly sampled ice data, new imputed environment, depth and doy

%LM_Model = 'Boost_RegAndImputedNoDoy'; LM_Model_title = 'Env, depth'; %regularly sampled ice data, imputed environment, depth no doy



Fram_0400_Init;

%Class model ID
R = load(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']));
Mdl1 = R.Mdl1;
bestParams = R.bestParams;
ClassifID = R.ClassifID;
%    'mask_seq2', 'cv_inner', 'cvErr' ,


predNames  = Mdl1.PredictorNames;

fGetPosterior = @(X) ml_returnPosterior(Mdl1, X);
classProbabil = fGetPosterior(tStPredictors(:, predNames));

%% Find mean chl-a
tChlTable = [tStPredictors(:, {'chlorophyll_a_0'}), tDMC(:, "ClustID")];

alpha = 0.25;
%tranf = @(X) X.^alpha; tranf_inv = @(X) X.^(1/alpha);
tranf = @(X) log10(X); tranf_inv = @(X) 10.^X;

tChlTable.chlorophyll_a_0 = tranf(tChlTable.chlorophyll_a_0);

quant = @(A) quantile(A,0.5);
%quant = @(A) mean(A, 'omitmissing');
tChlMean = grpstats(tChlTable(:, {'chlorophyll_a_0', 'ClustID'}), "ClustID", {'gname', quant}, 'VarNames',{'ClustID', 'gname', 'GroupCount', 'chlorophyll_a_0'});


%predict class
tChlTable.chlorophyll_a_predict =  classProbabil*tChlMean.chlorophyll_a_0;

lm = fitlm(tChlTable.chlorophyll_a_0, tChlTable.chlorophyll_a_predict,  'Intercept',true)


figure(1)
clf
s = scatter(tranf_inv(tChlTable.chlorophyll_a_0), tranf_inv(tChlTable.chlorophyll_a_predict),  "filled");
s.MarkerFaceAlpha = 0.3;


x_sorted = sort(tChlTable.chlorophyll_a_0)';
x_sorted = x_sorted(~isnan(x_sorted));
[y_fit, y_ci] = predict(lm, x_sorted', 'Alpha',0.01);
y_fit = tranf_inv(y_fit);
y_ci = tranf_inv(y_ci);
% Plot regression line
x_sorted_tranf = tranf_inv(x_sorted);
hold on
clr =  [0, 0, 0]; % [0.8000    0.2000         0];
if lm.Coefficients.pValue(2) <= 0.05
    plot(x_sorted_tranf, y_fit, '-', 'LineWidth', 3, 'Color',clr);
else
    plot(x_sorted_tranf, y_fit, '--', 'LineWidth', 1, 'Color',clr);
end
% Plot confidence interval
f = fill([x_sorted_tranf'; flipud(x_sorted_tranf')], [y_ci(:, 1); flipud(y_ci(:, 2))], clr, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
f.Annotation.LegendInformation.IconDisplayStyle ="off";

chl = linspace(min(tChlTable.chlorophyll_a_0), max(tChlTable.chlorophyll_a_0), 10);
%ln = plot(x_sorted_tranf, y_fit, 'Color','r', 'LineWidth',3);
plot(tranf_inv(chl), tranf_inv(chl), 'k--')
hold off
xlim([0.01, 4])
ylim([0.01, 4])
f_Set_scales({'logx', 'logy'})
xlabel('Observed')
ylabel('Predicted')
title('Chlorophyll a')

%%
% Number of bins along each axis
nBinsDoy = 10;
nBinsLong = 11;

% Short names for convenience
doy  = tStPredictors.Month_doy;
lon  = tStPredictors.Lon;
z    = tStPredictors.Depth;
chl  =  (tChlTable.chlorophyll_a_0);
chl_pred  = (tChlTable.chlorophyll_a_predict);
%chl_pred(isnan(chl)) = NaN;

% 1) Find min/max for axes
minDoy = 6;
maxDoy = 9.5;
minLon = -5;
maxLon = 10;

% 2) Construct bin edges (10x10 grid)
doyEdges = linspace(minDoy, maxDoy, nBinsDoy + 1);
lonEdges = linspace(minLon, maxLon, nBinsLong + 1);

% 3) Depth layers: [0–30), [30–60), [60–100)
layers = [  0  30;
    30  60;
    60 100];

% Preallocate array for mean chlorophyll:
% chlMean(:,:,k) –
chlMean = NaN(nBinsDoy, nBinsLong, size(layers,1));
chlPredMean = NaN(nBinsDoy, nBinsLong, size(layers,1));

%
for k = 1:size(layers,1)
    zmin = layers(k,1);
    zmax = layers(k,2);

    % Logical mask for current depth layer
    inLayer = (z >= zmin) & (z < zmax);

    doy_k = doy(inLayer);
    lon_k = lon(inLayer);
    chl_k = chl(inLayer);
    chl_pred_k = chl_pred(inLayer);

    % Skip if no data in this layer
    if isempty(doy_k)
        continue
    end

    % 4) Discretize DOY and Lon into bins
    %    doyBin, lonBin are bin indices from 1..nBins (or NaN if out of range)
    doyBin = discretize(doy_k, doyEdges);
    lonBin = discretize(lon_k, lonEdges);

    valid = ~isnan(doyBin) & ~isnan(lonBin);

    % 5) Compute mean chlorophyll in each (doy, lon) bin
    %    Empty bins will be filled with NaN (by default value)
    M = accumarray([doyBin(valid), lonBin(valid)], chl_k(valid), ...
        [nBinsDoy, nBinsLong], @median, NaN);

    chlMean(:,:,k) = M;

    M = accumarray([doyBin(valid), lonBin(valid)], chl_pred_k(valid), ...
        [nBinsDoy, nBinsLong], @median, NaN);
    chlPredMean(:,:,k) = M;
end
ca = [0.6, 2; ...
    0.2, 1;...
    0.1,1 ];

%
figure(2);
tiledlayout(3, 2, "TileSpacing", "compact", "Padding", "compact");
nLayers = size(chlMean, 3);
for k = 1:nLayers
    nexttile;
    [M, N] = size(chlMean(:,:,k));
    C_new = [chlMean(:,:,k) NaN(M, 1); NaN(1, N + 1)];
    pcolor(lonEdges, doyEdges, tranf_inv(C_new))
    shading flat
    colorbar
    caxis(ca(k, :));
    set(gca,'ColorScale','log')
    drawnow
    nexttile
    [M, N] = size(chlPredMean(:,:,k));
    C_new2 = [chlPredMean(:,:,k) NaN(M, 1); NaN(1, N + 1)];
    pcolor(lonEdges, doyEdges, tranf_inv(C_new2))
    shading flat
    caxis(ca(k, :))
    set(gca,'ColorScale','log')
    colorbar
end