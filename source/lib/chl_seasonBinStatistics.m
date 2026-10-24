function [stats, profiles] = chl_seasonBinStatistics(T, options)
%CHL_SEASONBINSTATISTICS Equal-month seasonal mean and profile percentiles.
% T must contain complete profiles, including observations below the layer.
% ProfileP05/ProfileP95 describe profile variability, not mean uncertainty.
arguments
    T table
    options.ProfileVar (1,:) char = 'ProfileID'
    options.MonthVar (1,:) char = 'Month_doy'
    options.DepthVar (1,:) char = 'Depth'
    options.ChlVar (1,:) char = 'chlorophyll_a_0'
    options.DepthRange (1,2) double = [0 30]
    options.DepthStep (1,1) double {mustBePositive} = 1
    options.MaxTopDepth (1,1) double = 10
    options.MinBottomDepth (1,1) double = 20
    options.SeasonMonths (1,:) double = 6:9
end

required = {options.ProfileVar, options.MonthVar, options.DepthVar, options.ChlVar};
assert(all(ismember(required, T.Properties.VariableNames)), ...
    'Missing required profile, month, depth, or chlorophyll column.');
assert(options.DepthRange(2) > options.DepthRange(1), 'DepthRange must increase.');
assert(options.MaxTopDepth >= options.DepthRange(1) && ...
    options.MaxTopDepth < options.MinBottomDepth && ...
    options.MinBottomDepth <= options.DepthRange(2), ...
    'Coverage thresholds must lie inside DepthRange and increase.');

stats = struct('MeanChl', NaN, 'ProfileP05', NaN, 'ProfileP95', NaN, ...
    'NProfiles', 0, 'NMonths', 0);
profiles = table(zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'ProfileGroup', 'Month', 'MeanChl'});
if isempty(T)
    return
end

groups = findgroups(T(:, options.ProfileVar));
depthGrid = unique([options.DepthRange(1):options.DepthStep:options.DepthRange(2), ...
    options.DepthRange(2)])';
for g = unique(groups(isfinite(groups)))'
    P = T(groups == g, :);
    months = unique(floor(P.(options.MonthVar)(isfinite(P.(options.MonthVar)))));
    if numel(months) ~= 1 || ~ismember(months, options.SeasonMonths)
        continue
    end
    valid = isfinite(P.(options.DepthVar)) & isfinite(P.(options.ChlVar));
    z = P.(options.DepthVar)(valid);
    c = P.(options.ChlVar)(valid);
    % Collapse duplicate depths before shape-preserving interpolation.
    [z, ~, depthGroup] = unique(z);
    if numel(z) < 2 || z(1) > options.MaxTopDepth || z(end) < options.MinBottomDepth
        continue
    end
    c = accumarray(depthGroup, c, [], @mean);
    % Clamp query depths to observations: constant extension at both ends.
    queryDepth = min(max(depthGrid, z(1)), z(end));
    cGrid = interp1(z, c, queryDepth, 'pchip');
    profileMean = trapz(depthGrid, cGrid) / diff(options.DepthRange);
    profiles(end + 1, :) = {g, months, profileMean};
end

if isempty(profiles)
    return
end
monthGroup = findgroups(profiles.Month);
monthlyMeans = splitapply(@mean, profiles.MeanChl, monthGroup);
stats.MeanChl = mean(monthlyMeans);
p = prctile(profiles.MeanChl, [5 95]);
stats.ProfileP05 = p(1);
stats.ProfileP95 = p(2);
stats.NProfiles = height(profiles);
stats.NMonths = numel(monthlyMeans);
end
