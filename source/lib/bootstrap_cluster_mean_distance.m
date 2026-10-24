function [meanDist, CI, bootMeans] = bootstrap_cluster_mean_distance(D, clusterID, k, nBoot)
% Bootstrap CI for mean within-cluster pairwise distance
%
% Inputs:
%   D         - NxN symmetric distance matrix
%   clusterID - Nx1 vector of cluster labels
%   k         - cluster label to analyze
%   nBoot     - number of bootstrap replicates
%
% Outputs:
%   meanDist  - observed mean within-cluster distance
%   CI        - 95% bootstrap confidence interval [low high]
%   bootMeans - bootstrap replicate means

    if nargin < 4
        nBoot = 1000;
    end

    % Indices of samples in cluster k
    idx = find(clusterID == k);
    n = numel(idx);

    if n < 2
        meanDist = NaN;
        CI = [NaN NaN];
        bootMeans = NaN(nBoot,1);
        return
    end

    % Observed mean distance in the original cluster
    Dk = D(idx, idx);
    maskOrig = triu(true(n), 1);
    dOrig = Dk(maskOrig);
    dOrig = dOrig(~isnan(dOrig));
    meanDist = mean(dOrig);

    % Bootstrap
    bootMeans = NaN(nBoot, 1);

    for b = 1:nBoot
        % Resample cluster members with replacement
        bootPos = randi(n, n, 1);   % positions within idx
        idxBoot = idx(bootPos);     % corresponding original sample indices

        % Bootstrap distance matrix
        Db = D(idxBoot, idxBoot);

        % Upper triangle
        mask = triu(true(n), 1);

        % Exclude artificial self-pairs caused by repeated sampling
        % i.e. pairs where the same original object was drawn twice
        sameObject = (idxBoot(:) == idxBoot(:)');
        mask = mask & ~sameObject;

        dBoot = Db(mask);
        dBoot = dBoot(~isnan(dBoot));

        if ~isempty(dBoot)
            bootMeans(b) = mean(dBoot);
        end
    end

    % Remove rare NaN replicates
    bootMeans = bootMeans(~isnan(bootMeans));

    % Percentile CI
    CI = prctile(bootMeans, [2.5 97.5]);
end