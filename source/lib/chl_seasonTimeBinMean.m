function [seasonalMean, nTimeBins, nMonths] = chl_seasonTimeBinMean(values, monthDoy, seasonStart, seasonEnd, monthBinStep)
%CHL_SEASONTIMEBINMEAN Average profile values, then nonempty time bins equally.
arguments
    values (:,1) double
    monthDoy (:,1) double
    seasonStart (1,1) double {mustBeFinite} = 6
    seasonEnd (1,1) double {mustBeFinite} = 10
    monthBinStep (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 0.5
end
assert(numel(values) == numel(monthDoy), 'Values and Month_doy must have equal lengths.');
assert(seasonEnd > seasonStart, 'Season end must exceed season start.');
edges = unique([seasonStart:monthBinStep:seasonEnd, seasonEnd]);
valid = isfinite(values) & isfinite(monthDoy) & ...
    monthDoy >= seasonStart & monthDoy < seasonEnd;
seasonalMean = NaN;
nTimeBins = 0;
nMonths = 0;
if ~any(valid)
    return
end
values = values(valid);
monthDoy = monthDoy(valid);
[~, ~, binID] = histcounts(monthDoy, edges);
group = findgroups(binID);
binMeans = splitapply(@mean, values, group);
seasonalMean = mean(binMeans);
nTimeBins = numel(binMeans);
nMonths = numel(unique(floor(monthDoy)));
end
