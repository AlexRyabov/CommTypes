%use
proj_Initialize;
rng(6)

%% Load data
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC.xlsx'));

%%
switch HGP.dmCluster.Method
    case 'rad' % Radius clustering
        clust_model = radius_cluster(tDMC{:, HGP.dmCluster.DMCoords}, 'Radius',HGP.dmCluster.rad_clust, 'Manner', 'approx', 'ReturnCenters', true);
        tDMC.ClustID = clust_model.labels;
    case 'hierarch' % Hierarchical clustering
        %load dm data
        tmp = load(fullfile(HGP.Data.Path2ResultsData, 'DM.mat'));
        DM = tmp.DM; clearvars tmp;
        Z = linkage(DM.simM.Dfull, 'complete');  % 'complete' or 'ward', ...
        labels = cluster(Z, 'maxclust', HGP.dmCluster.NumComponents);
        clust_model.gmm = 'hierarch';
        clust_model.K = HGP.dmCluster.NumComponents;
        clust_model.labels = labels;
        tDMC.ClustID = labels;
        [fg, MaxFigure] = addFindFigure('Dendrogram', MaxFigure, [-1, -1, 400 400], 'normal', 'on');
        dendrogram(Z)

    case 'gmm' %gaussian clustering
        X = tDMC{:, HGP.dmCluster.DMCoords};
        opts = statset('MaxIter', 1000, 'TolFun', 1e-6);
        gm = fitgmdist(X, HGP.dmCluster.NumComponents, ...
            'CovarianceType',  'diagonal', ...
            'SharedCovariance', true, ...
            'RegularizationValue', 1e-6, ...
            'Start','plus', ...
            'Replicates', 100, ...
            'Options', opts);
        clust_model.gmm = 'gmm';
        clust_model.K = gm.NumComponents;
        [P, ~] = posterior(gm, X);
        clust_model.gmmP = P;
        [~, labels] = max(P, [], 2);
        if  HGP.dmCluster.NumComponents == 9
            map = [7, 2, 8, 9, 6, 1, 3, 4, 5];
            labels_new = map(labels)';
        elseif HGP.dmCluster.NumComponents == 8
            map = [6, 4, 2, 5, 3, 1, 8, 7];
            labels_new = map(labels)';
            if strcmp(JobID, 'ASVNoMetazoaJBlocked')
                %remap types
                map2 = [7, 5, 4, 8, 2, 6, 3, 1];
                labels_new  = map2(labels_new)';
            end
        else
            labels_new = labels;
        end
        clearvars labels
        clust_model.labels = labels_new;
        clust_model.BIC = gm.BIC;
        clust_model.ICL = gm.BIC + sum(P(:) .* log(max(P(:), eps)));
        tDMC.ClustID = clust_model.labels;
    otherwise
        error('Underfined %s', HGP.dmCluster.Method)
end

writenewtable(tDMC, fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));

save(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'), 'clust_model');

%%
[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1, 1162 400], 'normal', 'on');
nexttile
clrs = f_Clrs_fresh(clust_model.K);
hold on
for k = 1:clust_model.K
    idx = (clust_model.labels == k);
    scatter3(tDMC{idx,1}, tDMC{idx,2}, tDMC{idx,3}, ...
        40, repmat(clrs(k,:), sum(idx),1), 'filled', 'MarkerEdgeColor', 'k');
end
legend('Location','best')
f_Lbls('DM_1', 'DM_2', 'DM_3')
hold off

nexttile
clrs = f_Clrs_fresh(clust_model.K);
hold on
for k = 1:clust_model.K
    idx = (clust_model.labels == k);
    scatter3(tDMC{idx,1}, tDMC{idx,2}, tDMC{idx,4}, ...
        40, repmat(clrs(k,:), sum(idx),1), 'filled', 'MarkerEdgeColor', 'k');
end
legend('Location','best')
hold off
f_Lbls('DM_1', 'DM_2', 'DM_4')

%% Jaccard distances and histograms: all pairs vs clusters

%% Pairwise distance distributions within community types

idxSamplCols = DI.idxSamplCols;
threshold    = 0.00;   % presence/absence threshold
useJaccard   = true;

if useJaccard
    % Binary matrix for Jaccard: rows = samples, columns = ASVs/features
    X = tSamplesNorm{:, idxSamplCols}';
    X = X > threshold;

    % All pairwise Jaccard distances: no diagonal, no duplicates
    dAll = pdist(X, 'jaccard');

    xLab  = 'Jaccard distance';
    edges = 0:0.01:1;
else
    % Distances in diffusion-map space
    X = tDMC{:, HGP.dmCluster.DMCoords};

    dAll = pdist(X, 'euclidean');

    xLab  = 'Diffusion-map distance';
    edges = linspace(min(dAll), max(dAll), 101);
end

dAll = dAll(~isnan(dAll));
mean_dist_all = mean(dAll, 'omitnan');

% Square distance matrix for cluster-based subsetting
D = squareform(dAll);
nObj = size(D, 1);

% Remove self-distances
D(1:nObj+1:end) = NaN;

% Cluster labels
clustIDs = sort(unique(tDMC.ClustID(:)))';
K = numel(clustIDs);

% Colours
if ~exist('clrs', 'var') || size(clrs, 1) < K
    clrs = f_Clrs_fresh(K);
end

% Output statistics
mean_dist_cluster = NaN(K, 1);
meanDist          = NaN(K, 1);
meanDistCI        = NaN(K, 2);
distReduct        = NaN(K, 1);
distReductCI      = NaN(K, 2);
OVL               = NaN(K, 1);

[fg, MaxFigure] = addFindFigure('JaccDistAllAndClusters', MaxFigure, ...
    [10, 20, 800, 800], 'normal', 'on');
clf(fg);
figure(fg);

nRows = ceil(K / 2);
tl = tiledlayout(fg, nRows, 2, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

for ii = 1:K

    k = clustIDs(ii);
    ind = (tDMC.ClustID == k);
    nK  = nnz(ind);

    ax = nexttile(tl);

    if nK < 2
        title(ax, sprintf('Type %d (no pairs)', k));
        axis(ax, 'off');
        continue;
    end

    % Distance submatrix within community type k
    Dk = D(ind, ind);

    % Upper triangle without diagonal = unique within-type pairs
    mask = triu(true(nK), 1);
    dK = Dk(mask);
    dK = dK(~isnan(dK));

    if isempty(dK)
        title(ax, sprintf('Type %d (no valid pairs)', k));
        axis(ax, 'off');
        continue;
    end

    % Mean within-type distance
    mean_dist_cluster(ii) = mean(dK, 'omitnan');

    % Bootstrap mean distance and reduction relative to all pairs
    [meanDist(ii), meanDistCI(ii, :), bootMeans] = ...
        bootstrap_cluster_mean_distance(D, tDMC.ClustID, k, 100);

    distRed_k = 1 - bootMeans / mean_dist_all;

    distReduct(ii)      = 1 - meanDist(ii) / mean_dist_all;
    distReductCI(ii, :) = prctile(distRed_k, [2.5 97.5]);

    % Plot full-dataset distance distribution
    histogram(ax, dAll, edges, ...
        'Normalization', 'pdf', ...
        'FaceColor', [1, 1, 1] * 0.8, ...
        'EdgeColor', 'none');

    hold(ax, 'on');

    % Plot within-type distance distribution
    histogram(ax, dK, edges, ...
        'Normalization', 'pdf', ...
        'FaceColor', clrs(ii, :), ...
        'FaceAlpha', 0.75, ...
        'EdgeColor', 'none');

    % Overlap coefficient as the shared area between two binned PDFs
    pdfAll = histcounts(dAll, edges, 'Normalization', 'pdf');
    pdfK   = histcounts(dK,   edges, 'Normalization', 'pdf');

    OVL(ii) = sum(min(pdfAll, pdfK) .* diff(edges));

    % Avoid tiny numerical deviations outside [0, 1]
    OVL(ii) = max(0, min(1, OVL(ii)));

    title(ax, sprintf('Type %d (overlap = %.0f%%)', ...
        k,  100 * OVL(ii)));

    xlabel(ax, xLab);
    ylabel(ax, 'PDF');

    if useJaccard
        xlim(ax, [0 1]);
        ylim(ax, [0 6.2]);
    end

    box(ax, 'on');
end

% Print summary statistics
fprintf('Mean distance across all pairs: %.3f\n', mean_dist_all);

summaryTable = table(clustIDs(:), mean_dist_cluster, meanDist, ...
    distReduct, distReductCI(:, 1), distReductCI(:, 2), OVL, ...
    'VariableNames', {'Type', 'MeanDistWithin', 'BootMeanDist', ...
                      'DistanceReduction', 'ReductionCI_low', ...
                      'ReductionCI_high', 'Overlap'});

disp(summaryTable);

max(mean_dist_cluster, [], 'omitnan')

%% plot silhouette 
[fg, MaxFigure] = addFindFigure('silhouette', MaxFigure, ...
    [10, 20, 800, 800], 'normal', 'on');

silhouette(tDMC{:, HGP.dmCluster.DMCoords}, tDMC.ClustID)

%% Number of samples of each community type by cruise
% DM rows and sample metadata follow the same sample order.
if height(tDMC) ~= height(tStPredictors)
    error('DM coordinates and sample metadata have different row counts.');
end

clusterIDs = sort(unique(tDMC.ClustID(:)))';
typeVars = cellstr("Type_" + string(clusterIDs));
tCruiseSamples = tStPredictors(:, {'ID', 'Cruise', 'Year', 'SamplingDate'});
tCruiseSamples.ClustID = categorical(tDMC.ClustID, clusterIDs, typeVars);

countsByCruiseType = grpstats(tCruiseSamples, ...
    {'Cruise', 'Year', 'ClustID'}, 'numel', 'DataVars', 'ID');
cruiseTypeCounts = unstack(countsByCruiseType, 'numel_ID', 'ClustID', ...
    'GroupingVariables', {'Cruise', 'Year'});

% A missing type in a cruise represents zero samples.
for iType = 1:numel(typeVars)
    if ~ismember(typeVars{iType}, cruiseTypeCounts.Properties.VariableNames)
        cruiseTypeCounts.(typeVars{iType}) = zeros(height(cruiseTypeCounts), 1);
    end
end
cruiseTypeCounts{:, typeVars} = fillmissing( ...
    cruiseTypeCounts{:, typeVars}, 'constant', 0);
cruiseTypeCounts = cruiseTypeCounts(:, [{'Cruise', 'Year'}, typeVars]);

firstSampleDates = grpstats(tCruiseSamples, {'Cruise', 'Year'}, ...
    'min', 'DataVars', 'SamplingDate');
firstSampleDates = firstSampleDates(:, ...
    {'Cruise', 'Year', 'min_SamplingDate'});
firstSampleDates.Properties.VariableNames{'min_SamplingDate'} = 'FirstSampleDate';

cruiseTypeCounts = join(firstSampleDates, cruiseTypeCounts, ...
    'Keys', {'Cruise', 'Year'});
cruiseTypeCounts = sortrows(cruiseTypeCounts, 'FirstSampleDate');
cruiseTypeCounts.Properties.RowNames = {};
disp(cruiseTypeCounts);
writenewtable(cruiseTypeCounts, fullfile(HGP.Data.Path2ResultsData, ...
    'Cruise_CommunityType_Counts.xlsx'));
