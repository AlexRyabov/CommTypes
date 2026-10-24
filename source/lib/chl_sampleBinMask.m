function [sampleRows, keepBins, counts] = chl_sampleBinMask(coordinates, edges, minSamples)
%CHL_SAMPLEBINMASK Exclude every sample in coordinate bins below the threshold.
arguments
    coordinates (:,1) double
    edges (1,:) double {mustBeFinite}
    minSamples (1,1) double {mustBeFinite, mustBePositive, mustBeInteger}
end
assert(numel(edges) >= 2 && all(diff(edges) > 0), 'Bin edges must increase.');
[~, ~, binID] = histcounts(coordinates, edges);
inRange = isfinite(coordinates) & binID > 0;
counts = accumarray(binID(inRange), ones(nnz(inRange), 1), [numel(edges)-1, 1]);
keepBins = counts >= minSamples;
sampleRows = false(size(coordinates));
sampleRows(inRange) = keepBins(binID(inRange));
end
