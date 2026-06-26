%use
%proj_Initialize;
rng(6)
%rng(42)

%% Do DM 
SpeciesAbundTransp = array2table(Data.dSamplesNorm');
spNames = cleanSpecNames(tSamplesNorm.TopName);
SpeciesAbundTransp.Properties.VariableNames = spNames;

%%find the best DM parameters for which we get the largest correlation with Depth and doy 
%DM = hg_cd_aux_FindBestDMParamters(k_min,Laplacian,sigma,alpha,DMvariables,SpeciesAbund,tStPredictors);
%metrics spearman, mcc, jaccard
%Laplacian='rownorm'
%k_min=55
k_min=60;
sigma=1.1;
alpha = 0.8;
Laplacian="rownorm";
%Laplacian="anisotropic";
Metric="jaccard";
%Metric="spearman";
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

tDMC = array2table(DM.DMC);
variableNames = arrayfun(@(x) sprintf('DM%d', x), 1:width(tDMC), 'UniformOutput', false);
tDMC.Properties.VariableNames = variableNames;



%% Define clusters 
HGP.dmCluster.NumComponents = 9;
%HGP.dmCluster.Method = 'gmm';
HGP.dmCluster.Method = 'medoids';
%HGP.dmCluster.Method = 'kmeans';
HGP.dmCluster.DMCoords = [1:6];
HGP.dmCluster.rad_clust = 0.3;
%%
clust_model.model = HGP.dmCluster.Method;
clust_model.BIC = NaN;
clust_model.ICL = NaN;

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
        clust_model.K = gm.NumComponents;
        [P, ~] = posterior(gm, X);
        clust_model.gmmP = P;
        [~, labels] = max(P, [], 2);
        if  HGP.dmCluster.NumComponents == 9
            map = [7, 2, 8, 9, 6, 1, 3, 4, 5];
            labels_new = map(labels)';
        elseif HGP.dmCluster.NumComponents == 8
            % order = [1, 2, 3, 4, 5, 6, 7, 8;
            % 4, 3, 7, 5, 6, 1, 8, 2]';
            %[~, ind] = sort(order(:, 2)); order(ind, :)
            map = [8, 5, 1, 3, 7, 6, 4, 2];
            %map = [6, 8, 3, 1, 4, 5, 2, 7];
            labels_new = map(labels)';
            % labels_new = labels;
        else
            labels_new = labels;
        end
        clearvars labels
        clust_model.labels = labels_new;
        clust_model.BIC = gm.BIC;
        clust_model.ICL = gm.BIC + sum(P(:) .* log(max(P(:), eps)));
        tDMC.ClustID = clust_model.labels;
   case 'kmeans' % k-means clustering
        X = tDMC{:, HGP.dmCluster.DMCoords};
        K = HGP.dmCluster.NumComponents;

        % Options for k-means
        opts = statset('MaxIter', 1000);

        % Hard clustering by k-means
        [labels, C, sumd] = kmeans(X, K, ...
            'Distance',   'sqeuclidean', ...
            'Replicates', 100, ...
            'Start',      'plus', ...
            'Options',    opts);
        clust_model.K = K;        
        clust_model.labels = labels;
        tDMC.ClustID = clust_model.labels;
    case 'medoids' % medoids
        idxSamplCols = DI.idxSamplCols;
        threshold    = 1e-5;   % presence/absence threshold

        % Binary matrix for Jaccard (rows = samples, columns = features)
        X = tSamplesNorm{:, idxSamplCols}';
        X = X > threshold;

        % All pairwise Jaccard distances (vector form)
        dAll = pdist(X, 'jaccard');           % [nPairs x 1]

        % Full distance matrix
        D = squareform(dAll);                 % nSamples x nSamples

        % Make sure diagonal is exactly zero
        n = size(D, 1);
        D(1:n+1:end) = 0;

        % Number of clusters (set as needed)
        K = HGP.dmCluster.NumComponents;      % or e.g. K = 8;

        % k-medoids clustering using precomputed distance matrix
        [labels, medoids, sumd] = kmedoids_precomputed(D, K, 100);
        %$kmedoids(D, K, ...
         %   'Distance', 'precomputed');

        % idx     - cluster index (1..K) for each sample
        % medoids - indices of medoid samples (rows of D / X)
        % sumd    - within-cluster sums of distances

        % Optionally save cluster IDs into your table
        tDMC.ClustID = labels;


     
        clust_model.K = K;        
        clust_model.labels = labels;
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



%%
%% Jaccard distances and histograms: all pairs vs clusters

idxSamplCols = DI.idxSamplCols;
threshold    = 1e-5;   % presence/absence threshold

% Binary matrix for Jaccard (rows = samples, columns = features)
X = tSamplesNorm{:, idxSamplCols}';
X = X > threshold;

% All pairwise Jaccard distances (vector, no diagonal, no duplicates)
dAll = pdist(X, 'jaccard');           % [nPairs x 1]


mean_dist_all = mean(dAll);

% Square distance matrix for cluster-based subsetting
D = squareform(dAll);                 % nObj x nObj
nObj = size(D, 1);

% Remove self-distances from further calculations
D(1:nObj+1:end) = NaN;

K = length(unique(tDMC.ClustID));
mean_dist_cluster = NaN(K, 1);

edges = 0:0.01:1;

[fg, MaxFigure] = addFindFigure('JaccDistAllAndClusters', MaxFigure, ...
    [10, 20, 800, 800], 'normal', 'on');
clf(fg);
figure(fg);

% Layout for cluster panels
nRows = ceil(K/2);
tiledlayout(fg, nRows, 2, 'Padding','compact', 'TileSpacing','compact');

for k = 1:K
    ind = (tDMC.ClustID == k);  % logical index of samples in cluster k
    nK  = nnz(ind);

    nexttile;

    if nK < 2
        % Not enough samples in this cluster to form pairs
        title(sprintf('Cluster %d (n = %d, no pairs)', k, nK));
        axis off;
        continue;
    end

    % Distance submatrix within cluster k
    Dk = D(ind, ind);           % nK x nK

    % Upper triangle without diagonal = unique pairs
    mask = triu(true(nK), 1);
    dK = Dk(mask);

    % Mean distance within cluster k
    mean_dist_cluster(k) = mean(dK, 'omitnan');

    % Background: all-pairs distribution (grey)
    histogram(dAll, edges, ...
        'Normalization','pdf', ...
        'FaceColor',[1, 1, 1]*0.8, ...
        'EdgeColor','none');
    hold on;

    % Overlay: within-cluster pairs (colored)
    histogram(dK, edges, ...
        'Normalization','pdf', ...
        'FaceColor',clrs(k,:), ...
        'EdgeColor','none');

    title(sprintf('Cluster %d (n = %d, mean = %.3f)', ...
        k, nK, mean_dist_cluster(k)));
    xlabel('Jaccard distance');
    ylabel('PDF');
    box on;
end

% Print summary statistics to command window
fprintf('Mean Jaccard distance (all pairs): %.3f\n', mean_dist_all);
disp(table((1:K).', mean_dist_cluster, ...
    'VariableNames', {'Cluster', 'MeanDistWithin'}));
max(mean_dist_cluster)

%%
