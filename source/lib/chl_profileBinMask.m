function [profileRows, sampleRows, keepBins, profileCounts, sampleCounts] = ...
    chl_profileBinMask(profileCoordinates, sampleCoordinates, edges, minProfiles)
%CHL_PROFILEBINMASK Keep bins by accepted-profile count, not sample count.
arguments
    profileCoordinates (:,1) double
    sampleCoordinates (:,1) double
    edges (1,:) double {mustBeFinite}
    minProfiles (1,1) double {mustBeFinite, mustBePositive, mustBeInteger}
end
[profileRows, keepBins, profileCounts] = chl_sampleBinMask( ...
    profileCoordinates, edges, minProfiles);
[~, ~, sampleCounts] = chl_sampleBinMask(sampleCoordinates, edges, 1);
[~, ~, sampleBin] = histcounts(sampleCoordinates, edges);
sampleRows = isfinite(sampleCoordinates) & sampleBin > 0;
sampleRows(sampleRows) = keepBins(sampleBin(sampleRows));
end
