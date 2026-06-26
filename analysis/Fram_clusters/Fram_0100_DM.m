%% Diffusion Map of Hausgarten samples with a new DM 
% Alex Ryabov, 2025
%We want to plot DM of samples, cluster the samples and find typical
%environmental conditions to characterize each cluster. 
%further, each cluster can be characterized with its specific community
%composition at the genera/class/order levels 
%the idea it to find the most important drivers of the communities 
%Plan
%1. DM
%2a Clusters of DM
%2b Characterize typical communities within each cluster. Which species
%dominate, what is the diversity +- within the cluster. 
%2c Try to find common species across the clusters and remove them, to do
%it find average species composition for each cluster. 
%3. Random forest or something similar to predict the cluster based on the
%environment 
%4. Partial dependencies, of comm

%Clean variables, close windows, initialize all project constants and parameters
%clearvars
proj_Initialize;

%% Load data 


%%
grpstats(tStPredictors(:, {'Year', 'MonthNum'}), {'Year', 'MonthNum'}, ["numel" ])

% Find the histogram of average species abundance 
SpeciesAbundTransp = array2table(Data.dSamplesNorm');
spNames = cleanSpecNames(tSamplesNorm.TopName);
SpeciesAbundTransp.Properties.VariableNames = spNames;
hg_cd_aux_plotspeceis(SpeciesAbundTransp, 'RankAbund2', 101);
%f_Lbls('')
%% DM 

%%find the best DM parameters for which we get the largest correlation with Depth and doy 
%DM = hg_cd_aux_FindBestDMParamters(k_min,Laplacian,sigma,alpha,DMvariables,SpeciesAbund,tStPredictors);
%metrics spearman, mcc, jaccard
%Laplacian='rownorm'
%k_min=60;
k_min=60;
sigma=1.3;
alpha = 1;
Laplacian="rownorm";
%Metric="jaccard";
Metric=HGP.dm.SimilarityMetric;
%Metric="SERa";
%Metric="dist_BrayCurtis";
%Metric="morisitaHorn";
%Metric="mcc";
DMvariables=15;

mask = [];
%mask = diffmap_mask(tStPredictors(:, ["Year", "HGStation"]));
%mask = diffmap_mask(tStPredictors(:, ["Year"]));

options = diffmap_options('k_min',k_min, 'Laplacian',Laplacian,  'Metric',Metric, 'sigma',sigma, 'alpha',alpha, 'DMvariables',DMvariables, 'mask', mask);
DM = diffmap(SpeciesAbundTransp{:, :}, options);
DM.DMC(:, 1) = -DM.DMC(:, 1);
DM.DMC(:, 2) = -DM.DMC(:, 2);
hg_dm_plotDM1(DM,tStPredictors);
disp('All species')
%% find and plot correlations
drCols = ["doy", "mod_aq_par", "temp", "mod_aq_sst", "sea_ice_monthly",  "sea_ice_now",  "salin", "mld", "Lon" , "Lat",  "chlorophyll_a_0", "mod_aq_par_depth", "phosphate", "nitrate",   "silicate", "DIN", "Depth", "DIN_100_200"];
dmCols = ["DMA1", "DMA2", "DMA3", "DMA4", "DMA5", "DMA6", "DMA7", "DMA8", "DMA9", "DMA10"];
maxDMaxes = 6;

drColsPlot = regexprep(drCols,'_',' ','emptymatch');
dmColsPlot = regexprep(dmCols,'_',' ','emptymatch');

[rho, p] = corr(DM.DMC(:, :), tStPredictors{:, drCols}, "Type","Spearman", "Rows","pairwise");
rho(p>0.01) = NaN;
[fg, MaxFigure] = addFindFigure('DiffMapCorr', MaxFigure, [10, 20, 1305 565], 'normal', 'on');
clf
t = tiledlayout(1, 3)
nexttile(1)
hm = heatmap(dmColsPlot(1:maxDMaxes), drColsPlot , rho(1:maxDMaxes, :)');
hm.MissingDataColor = [1, 1, 1]*0.9;
hm.CellLabelFormat = '%0.2g';
colormap(redbluecmap); caxis([-1, 1]);


%%Find correlations between predictors 
[rhoEnv, pEnv] = corr( tStPredictors{:, drCols}, tStPredictors{:, drCols}, "Type","Spearman", "Rows","pairwise");

%Sort them by hierarchical clustering 
D = 1 - abs(rhoEnv);                  % distance matrix
D(pEnv>0.0001) = 10;

Z = linkage(squareform(D),"complete"); % hierarchical clustering
%ord = optimalleaforder(Z, squareform(D)); % better visual ordering

%set nonsignificat elements to zero 
rhoEnv(pEnv>0.0001) = NaN;

nexttile([1, 2])
hm = heatmap(drColsPlot(ord), drColsPlot(ord), rhoEnv(ord, ord));
hm.MissingDataColor = [1, 1, 1]*0.9;
hm.CellLabelFormat = '%0.2g';
colormap(redbluecmap); caxis([-1, 1]);
AddLetters2Plots( 'HShift', [-0.08])

%%
% [rho, p] = corr(DM.DMC(:, :), tStPredictors{:, 16:end}, ...
%     "Type","Spearman", "Rows","pairwise");
% rho(p>0.000001) = NaN;
% 
[bcR, p, T, df] = bcdistcorr(DM.DMC(:, 1:5), tStPredictors{:, ["doy", "Depth"]});

% 
tDMC = array2table(DM.DMC);
variableNames = arrayfun(@(x) sprintf('DM%d', x), 1:width(tDMC), 'UniformOutput', false);
tDMC.Properties.VariableNames = variableNames;
mkdir(HGP.Data.Path2ResultsData);
writenewtable(tDMC, fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC.xlsx'));
save(fullfile(HGP.Data.Path2ResultsData, 'DM.mat'), "DM");


%% plot Metric distance vs diffusion distance 

metricD = pdist(SpeciesAbundTransp{:, :}>0, Metric);           % [nPairs x 1]
metricD = squareform(metricD);                 % nObj x nObj

%metricD = 1-DM.simM.S;

DM_D = pdist(DM.DMC, 'euclidean'); 
DM_D = squareform(DM_D); 

[fg, MaxFigure] = addFindFigure('DM vs metric2', MaxFigure, ...
    [10, 20, 800, 800], 'normal', 'on');

ind = metricD(:) < 1;
metricD = metricD(ind );
DM_D = DM_D(ind );
s = scatter(metricD(:), DM_D(:), 'filled')
s.MarkerFaceAlpha = 0.1;
fitlm(metricD(:), DM_D(:))
xlabel('Metric distance');
ylabel('DM distance');