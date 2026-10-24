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
proj_Initialize;

%% Load data 


%%
grpstats(tStPredictors(:, {'Year', 'MonthNum'}), {'Year', 'MonthNum'}, ["numel" ])

% Find the histogram of average species abundance 
SpeciesAbundTransp = array2table(Data.dSamplesNorm');
spNames = cleanSpecNames(tSamplesNorm.TopName);
SpeciesAbundTransp.Properties.VariableNames = spNames;
hg_cd_aux_plotspeceis(SpeciesAbundTransp, 'RankAbund2', 101);
%% DM 

%%find the best DM parameters for which we get the largest correlation with Depth and doy 
k_min=61;
sigma=1.3;
alpha = 1;
Laplacian="rownorm";
Metric="jaccard";
DMvariables=15;

if HGP.dm.BlockByYears
    mask = diffmap_mask(tStPredictors(:, ["Year"]));
else 
    mask = [];
end

options = diffmap_options('k_min',k_min, 'Laplacian',Laplacian, ...
 'Metric',Metric, 'sigma',sigma, 'alpha',alpha, 'DMvariables',DMvariables, ...
 'mask', mask, 'diffDistance', true);
DM = diffmap(SpeciesAbundTransp{:, :}, options);
DM.DMC(:, 1) = -DM.DMC(:, 1);
DM.DMC(:, 2) = DM.DMC(:, 2);
hg_dm_plotDM1(DM,tStPredictors);
disp('All species')

%% find and plot correlations
drCols = ["doy", "mod_aq_par", "temp", "mod_aq_sst", "sea_ice_monthly",  "sea_ice_now",  "salin", "mld", "Lon" , "Lat",  "chlorophyll_a_0", "mod_aq_par_depth", "phosphate", "nitrate",   "silicate", "DIN", "Depth", "DIN_100_200"];
dmCols = ["DMA1", "DMA2", "DMA3", "DMA4", "DMA5", "DMA6", "DMA7", "DMA8", "DMA9", "DMA10"];
maxDMaxes = 6;
heatmapFontSize = 10;

drColsPlot = renamePredictors(drCols);
% Diffusion-map labels
dmColsPlot = "DM" + string(1:numel(dmCols));


[rho, p] = corr(DM.DMC(:, :), tStPredictors{:, drCols}, "Type","Spearman", "Rows","pairwise");
rho(p>0.0001) = NaN;


%%Find correlations between predictors 
[rhoEnv, pEnv] = corr( tStPredictors{:, drCols}, tStPredictors{:, drCols}, "Type","Spearman", "Rows","pairwise");

%Sort them by hierarchical clustering 
D = 1 - abs(rhoEnv);                  % distance matrix
D(pEnv>0.0001) = 10;

Z = linkage(squareform(D),"complete"); % hierarchical clustering
ord = optimalleaforder(Z, squareform(D)); % better visual ordering
predOrd = ord(:)';
predLabels = drColsPlot(predOrd);

%set nonsignificat elements to zero 
rhoEnv(pEnv>0.0001) = NaN;

[fg, MaxFigure] = addFindFigure('DiffMapCorr', MaxFigure, [10, 20, 1658 573], 'normal', 'on');
clf
t = tiledlayout(1, 4);
t.TileSpacing = 'tight';
t.Padding = 'compact';

nexttile(1);
hm = heatmap(dmColsPlot(1:maxDMaxes), predLabels, rho(1:maxDMaxes, predOrd)');
hm.FontSize = heatmapFontSize;
hm.MissingDataColor = [1, 1, 1]*0.9;
hm.CellLabelFormat = '%0.2g';
colormap(redbluecmap); caxis([-1, 1]);

nexttile(2, [1, 2]);
hm = heatmap(predLabels, predLabels, rhoEnv(predOrd, predOrd));
hm.FontSize = heatmapFontSize;
hm.MissingDataColor = [1, 1, 1]*0.9;
hm.CellLabelFormat = '%0.2g';
colormap(redbluecmap); caxis([-1, 1]);

nexttile(4);
ZPlot = Z;
ZPlot(:, 3) = log10(1 + ZPlot(:, 3));
dendrogram(ZPlot, 0, 'Reorder', predOrd, 'Orientation', 'right');
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', [], 'YTickLabel', []);
ylim([0.5, length(predLabels)+0.5])
xlabel('log_{10}(1 + distance)');
AddLetters2Plots(fg, {'a', 'b', 'c'}, 'HShift', -0.08);

%%
[bcR, p, T, df] = bcdistcorr(DM.DMC(:, 1:5), tStPredictors{:, ["doy", "Depth"]});

tDMC = array2table(DM.DMC);
variableNames = arrayfun(@(x) sprintf('DM%d', x), 1:width(tDMC), 'UniformOutput', false);
tDMC.Properties.VariableNames = variableNames;
mkdir(HGP.Data.Path2ResultsData);
writenewtable(tDMC, fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC.xlsx'));
save(fullfile(HGP.Data.Path2ResultsData, 'DM.mat'), "DM");


%% plot Metric distance vs diffusion distance 

metricD = pdist(SpeciesAbundTransp{:, :}>0, Metric);           % [nPairs x 1]
metricD = squareform(metricD);                 % nObj x nObj


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


%% Find fractal dimension 
[fg, MaxFigure] = addFindFigure('Correlation dimension estimate', MaxFigure, [10, 20, 378 290], 'normal', 'on');
[d, fitmodel] = f_CorrelationDimension(DM.Ddiff  , [0.1, 1])
xlim([0.08, 1.2])
ax = gca;
ax.Title.String = '';
