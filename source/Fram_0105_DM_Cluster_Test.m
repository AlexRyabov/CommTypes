%% Test of diffusion mapping and clustering with DM fully independent of environmental variables.
% Environmental variables are used only to define the year mask and to
% evaluate the finished map; they are not used as diffusion-map features.

proj_Initialize;

%% Diffusion map
SpeciesAbundTransp = array2table(Data.dSamplesNorm');
spNames = cleanSpecNames(tSamplesNorm.TopName);
SpeciesAbundTransp.Properties.VariableNames = spNames;

k_min = 28;
sigma = 1.3;
alpha = 1;
Laplacian = "rownorm";
Metric = "jaccard";
DMvariables = 15;
mask = diffmap_mask(tStPredictors(:, ["Year"]));

options = diffmap_options('k_min', k_min, 'Laplacian', Laplacian, ...
    'Metric', Metric, 'sigma', sigma, 'alpha', alpha, ...
    'DMvariables', DMvariables, 'mask', mask, 'diffDistance', true);
DM = diffmap(SpeciesAbundTransp{:, :}, options);
DM.DMC(:, 1) = -DM.DMC(:, 1);

%% Correlations of DM axes with environmental conditions
drCols = ["doy", "mod_aq_par", "temp", "mod_aq_sst", ...
    "sea_ice_monthly", "sea_ice_now", "salin", "mld", "Lon", "Lat", ...
    "chlorophyll_a_0", "mod_aq_par_depth", "phosphate", "nitrate", ...
    "silicate", "DIN", "Depth", "DIN_100_200"];
maxDMaxes = 6;
[rho, p] = corr(DM.DMC, tStPredictors{:, drCols}, ...
    'Type', 'Spearman', 'Rows', 'pairwise');
rho(p > 0.0001) = NaN;

% Use the same environmental-variable order as in Fram_0100_DM.
[rhoEnv, pEnv] = corr(tStPredictors{:, drCols}, ...
    tStPredictors{:, drCols}, 'Type', 'Spearman', 'Rows', 'pairwise');
distEnv = 1 - abs(rhoEnv);
distEnv(pEnv > 0.0001) = 10;
Z = linkage(squareform(distEnv), 'complete');
predOrd = optimalleaforder(Z, squareform(distEnv));
predLabels = renamePredictors(drCols);
predLabels = predLabels(predOrd);

%% Correlation dimension of the diffusion map
% Same correlation-integral fit and radius interval as Fram_0100_DM,
% without opening an additional figure.
D_upper = triu(DM.Ddiff, 1);
dists = nonzeros(D_upper);
r_values = logspace(log10(min(dists)), log10(max(dists)), 100);
C = arrayfun(@(r) sum(dists < r) / (size(DM.Ddiff, 1) * ...
    (size(DM.Ddiff, 1) - 1) / 2), r_values);
idx2fit = r_values > 0.1 & r_values < 1 & C > 0;
fitmodel = fitlm(log10(r_values(idx2fit))', log10(C(idx2fit))');
d = fitmodel.Coefficients.Estimate(2);
se = fitmodel.Coefficients.SE(2);
disp("Correlation dimension: " + string(d) + " +/- " + string(se));

%% Eight GMM clusters, using the settings from Fram_0250_Cluster
rng(6);
K = HGP.dmCluster.NumComponents;
assert(strcmp(HGP.dmCluster.Method, 'gmm') && K == 8, ...
    'This test expects the default eight-component GMM clustering.');
X = DM.DMC(:, HGP.dmCluster.DMCoords);
opts = statset('MaxIter', 1000, 'TolFun', 1e-6);
gm = fitgmdist(X, K, 'CovarianceType', 'diagonal', ...
    'SharedCovariance', true, 'RegularizationValue', 1e-6, ...
    'Start', 'plus', 'Replicates', 100, 'Options', opts);
[P, ~] = posterior(gm, X);
[~, rawLabels] = max(P, [], 2);

% Match each new cluster to one default cluster using correlations of
% sample-membership indicators. The one-to-one assignment maximizes the
% total correlation, so matched IDs use the same colors as the default map.
defaultFile = fullfile(HGP.Data.Path2ResultsData, ...
    'Samples_DMC_clustID.xlsx');
assert(isfile(defaultFile), ...
    'Run Fram_0250_Cluster on the default DM before this test.');
defaultDMC = readtable(defaultFile);
assert(height(defaultDMC) == numel(rawLabels), ...
    'The default clusters and new DM have different sample counts.');
defaultLabels = defaultDMC.ClustID;
assert(all(ismember(defaultLabels, 1:K)), ...
    'The default cluster file must contain labels 1 through 8.');

newMembership = double(rawLabels == 1:K);
defaultMembership = double(defaultLabels == 1:K);
clusterCorr = corr(newMembership, defaultMembership);
assert(all(isfinite(clusterCorr), 'all'), ...
    'All eight new and default clusters must contain samples.');
assignments = perms(1:K);
rowIndices = repmat(1:K, size(assignments, 1), 1);
scores = sum(clusterCorr(sub2ind([K, K], rowIndices, assignments)), 2);
[~, bestAssignment] = max(scores);
clusterMap = assignments(bestAssignment, :);
ClustID = reshape(clusterMap(rawLabels), [], 1);

overlap = accumarray([rawLabels, defaultLabels], 1, [K K]);
disp(table((1:K)', clusterMap', diag(overlap(:, clusterMap)), ...
    'VariableNames', {'NewCluster', 'DefaultCluster', 'SharedSamples'}));
disp('New-cluster rows by default-cluster columns (sample overlap):');
disp(array2table(overlap, 'VariableNames', "Default" + string(1:K), ...
    'RowNames', cellstr("New" + string(1:K))));

%% Three-panel comparison

[fg, MaxFigure] = addFindFigure('DM cluster test', MaxFigure, [50 50  1357 387], 'normal', 'on');


tl = tiledlayout(fg, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
mark_size = 30;
az_el = [0 90];

ax1 = nexttile(tl, 1);
scatter3(ax1, DM.DMC(:, 1), DM.DMC(:, 2), DM.DMC(:, 3), ...
    mark_size, DM.DMC(:, 3), 'filled', 'MarkerEdgeColor', 'k');
colormap(ax1, parula);
colorbar(ax1);
view(ax1, az_el);
xlabel(ax1, 'DM_1'); ylabel(ax1, 'DM_2'); zlabel(ax1, 'DM_3');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
clrs = f_Clrs_fresh(K);
for k = 1:K
    idx = ClustID == k;
    s = scatter3(ax2, DM.DMC(idx, 1), DM.DMC(idx, 2), ...
        DM.DMC(idx, 3), mark_size, clrs(k, :), 'filled', ...
        'MarkerEdgeColor', 'k');
    s.DisplayName = sprintf('C %d', k);
end
hold(ax2, 'off');
view(ax2, az_el);
xlabel(ax2, 'DM_1'); ylabel(ax2, 'DM_2'); zlabel(ax2, 'DM_3');
legend(ax2, 'Location', 'bestoutside');

nexttile(tl, 3);
hm = heatmap("DM" + string(1:maxDMaxes), predLabels, ...
    rho(1:maxDMaxes, predOrd)');
hm.FontSize = 10;
hm.MissingDataColor = [1 1 1] * 0.9;
hm.CellLabelFormat = '%0.2g';
hm.Colormap = redbluecmap;
hm.ColorLimits = [-1 1];
AddLetters2Plots(fg, {'a', 'b', 'c'}, HShift=-0.04, VShift=-0.01)

%% Pairwise Jaccard distances within each community type
% Adapted from the corresponding section of Fram_0250_Cluster.
threshold = 0;
Xbinary = tSamplesNorm{:, DI.idxSamplCols}' > threshold;
dAllPairs = pdist(Xbinary, 'jaccard');
distMatrix = squareform(dAllPairs);
dAll = dAllPairs(~isnan(dAllPairs));
mean_dist_all = mean(dAll, 'omitnan');
nObj = size(distMatrix, 1);
distMatrix(1:nObj+1:end) = NaN;
edges = 0:0.01:1;
pdfAll = histcounts(dAll, edges, 'Normalization', 'pdf');

mean_dist_cluster = NaN(K, 1);
meanDist = NaN(K, 1);
distReduct = NaN(K, 1);
distReductCI = NaN(K, 2);
OVL = NaN(K, 1);


[fgDist, MaxFigure] = addFindFigure('JaccDistAllAndClusters', MaxFigure, ...
    [10, 20, 800, 800], 'normal', 'on');
clf(fgDist);
figure(fgDist);

nRows = ceil(K / 2);
tlDist = tiledlayout(fgDist, nRows, 2, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

for k = 1:K
    idx = ClustID == k;
    nK = nnz(idx);
    ax = nexttile(tlDist);
    if nK < 2
        title(ax, sprintf('Type %d (no pairs)', k));
        axis(ax, 'off');
        continue;
    end

    Dk = distMatrix(idx, idx);
    dK = Dk(triu(true(nK), 1));
    dK = dK(~isnan(dK));
    if isempty(dK)
        title(ax, sprintf('Type %d (no valid pairs)', k));
        axis(ax, 'off');
        continue;
    end

    mean_dist_cluster(k) = mean(dK, 'omitnan');
    [meanDist(k), ~, bootMeans] = ...
        bootstrap_cluster_mean_distance(distMatrix, ClustID, k, 100);
    distRed_k = 1 - bootMeans / mean_dist_all;
    distReduct(k) = 1 - meanDist(k) / mean_dist_all;
    distReductCI(k, :) = prctile(distRed_k, [2.5 97.5]);

    histogram(ax, dAll, edges, 'Normalization', 'pdf', ...
        'FaceColor', [1 1 1] * 0.8, 'EdgeColor', 'none');
    hold(ax, 'on');
    histogram(ax, dK, edges, 'Normalization', 'pdf', ...
        'FaceColor', clrs(k, :), 'FaceAlpha', 0.75, 'EdgeColor', 'none');
    hold(ax, 'off');

    pdfK = histcounts(dK, edges, 'Normalization', 'pdf');
    OVL(k) = sum(min(pdfAll, pdfK) .* diff(edges));
    OVL(k) = max(0, min(1, OVL(k)));
    title(ax, sprintf('Type %d (overlap = %.0f%%)', ...
        k, 100 * OVL(k)));
    xlabel(ax, 'Jaccard distance');
    ylabel(ax, 'PDF');
    xlim(ax, [0 1]);
    ylim(ax, [0 6.2]);
    box(ax, 'on');
end

fprintf('Mean distance across all pairs: %.3f\n', mean_dist_all);
summaryTable = table((1:K)', mean_dist_cluster, meanDist, ...
    distReduct, distReductCI(:, 1), distReductCI(:, 2), OVL, ...
    'VariableNames', {'Type', 'MeanDistWithin', 'BootMeanDist', ...
    'DistanceReduction', 'ReductionCI_low', 'ReductionCI_high', 'Overlap'});
disp(summaryTable);
