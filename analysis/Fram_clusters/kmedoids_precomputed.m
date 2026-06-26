function [labels, medoids, sumd] = kmedoids_precomputed(D, K, maxIter)
%KMEDOIDS_PRECOMPUTED Simple k-medoids using a precomputed distance matrix.
%
%   [labels, medoids, sumd] = kmedoids_precomputed(D, K, maxIter)
%
%   D       - n x n distance matrix, D(i,i) can be zero (self-distance).
%   K       - number of clusters.
%   maxIter - maximum number of iterations.
%
%   labels  - n x 1 vector, cluster index for each object (1..K).
%   medoids - K x 1 vector, indices of medoid objects (rows of D).
%   sumd    - K x 1 vector, within-cluster sum of distances.

    if nargin < 3
        maxIter = 100;
    end

    % Basic checks
    if size(D,1) ~= size(D,2)
        error('D must be a square distance matrix.');
    end

    n = size(D, 1);

    % Initialize medoids randomly without replacement
    medoids = randperm(n, K)';

    labels = zeros(n,1);

    for iter = 1:maxIter
        % Assignment step: assign each point to the nearest medoid
        for i = 1:n
            [~, labels(i)] = min(D(i, medoids));
        end

        oldMedoids = medoids;

        % Update step: for each cluster choose the best medoid
        for k = 1:K
            clusterInd = find(labels == k);

            % If cluster is empty, reinitialize its medoid randomly
            if isempty(clusterInd)
                medoids(k) = randi(n);
                continue;
            end

            % Distances within the cluster
            Dk = D(clusterInd, clusterInd);  % size: nk x nk

            % For each candidate medoid in the cluster compute total distance
            totalDist = sum(Dk, 2);          % nk x 1

            % Choose the point with minimal total distance as new medoid
            [~, minIdxLocal] = min(totalDist);
            medoids(k) = clusterInd(minIdxLocal);
        end

        % Check for convergence: if medoids do not change, stop
        if isequal(medoids, oldMedoids)
            break;
        end
    end

    % Final assignment with converged medoids
    for i = 1:n
        [~, labels(i)] = min(D(i, medoids));
    end

    % Within-cluster sums of distances
    sumd = zeros(K,1);
    for k = 1:K
        clusterInd = find(labels == k);
        if isempty(clusterInd)
            sumd(k) = NaN;
        else
            mk = medoids(k);
            sumd(k) = sum(D(clusterInd, mk));
        end
    end
end
