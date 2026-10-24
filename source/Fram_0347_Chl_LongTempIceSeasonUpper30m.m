%% Seasonal upper-30 m chlorophyll and sea ice by longitude and temperature
% Integrate pchip chlorophyll profiles over 0-30 m. Average profiles within
% Month_doy intervals, then average nonempty intervals with equal weight.
% Configurable profile percentiles describe variability, not a confidence interval.

proj_Initialize;

if ~exist('Data', 'var') || ~isfield(Data, 'tStPredictors')
    [Data, ~] = dt_loadData(HGP, []);
end
tStPredictors = Data.tStPredictors;

keepRows = tStPredictors.chlorophyll_a_0 < 1117.4;
tStPredictors = tStPredictors(keepRows, :);

%% Settings
chlVar = 'chlorophyll_a_0';
iceVar = 'sea_ice_now';
monthVar = 'Month_doy';
depthVar = 'Depth';
profileVars = {'Cruise', 'Year', 'doy', 'HGStation'};
seasonStart = 6;
seasonEnd = 9;
monthBinStep = 1/6; % any positive averaging interval in Month_doy units
validateattributes(monthBinStep, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
depthRange = [0 30];
depthStep = 1;
maxTopDepth = 10;
minBottomDepth = 20;
lonBinWidth = 1.5;
tempBinWidth = 0.75;
plotCommunityBubbles = true; % true shows chlorophyll bubbles by community type
bubbleAreaPerSample = 180 / 50; % bubble area in square points per sample
minProfilesPerBin = 4; % exclude bins with fewer accepted chlorophyll profiles
profilePercentiles = [20 80]; % lower/upper envelope percentiles, in percent
validateattributes(profilePercentiles, {'numeric'}, ...
    {'real', 'finite', 'size', [1 2], '>=', 0, '<=', 100});
assert(profilePercentiles(1) < profilePercentiles(2), ...
    'profilePercentiles must contain increasing lower and upper bounds.');

if plotCommunityBubbles
    tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
    tDMC = tDMC(keepRows, :);
end

requiredVars = [{chlVar, iceVar, monthVar, depthVar, 'Lon', 'temp'}, profileVars];
missingVars = requiredVars(~ismember(requiredVars, tStPredictors.Properties.VariableNames));
if ~isempty(missingVars)
    error('Missing required tStPredictors variable(s): %s', strjoin(missingVars, ', '));
end
tStPredictors.ProfileID = findgroups(tStPredictors(:, profileVars));
seasonRows = tStPredictors.(monthVar) >= seasonStart & ...
    tStPredictors.(monthVar) < seasonEnd & ...
    isfinite(tStPredictors.ProfileID) & ...
    isfinite(tStPredictors.(monthVar)) & ...
    isfinite(tStPredictors.(depthVar)) & ...
    isfinite(tStPredictors.(chlVar));
layerRows = seasonRows & tStPredictors.(depthVar) >= depthRange(1) & ...
    tStPredictors.(depthVar) <= depthRange(2);
T = tStPredictors(seasonRows, :);

%% Use the same accepted chlorophyll profiles for both variables
[~, integratedProfiles] = chl_seasonBinStatistics(T, ...
    MonthVar=monthVar, DepthVar=depthVar, ChlVar=chlVar, ...
    DepthRange=depthRange, DepthStep=depthStep, ...
    MaxTopDepth=maxTopDepth, MinBottomDepth=minBottomDepth, ...
    SeasonMonths=seasonStart:(seasonEnd - 1));
if isempty(integratedProfiles)
    error('No seasonal chlorophyll profiles meet the depth coverage requirements.');
end
[~, profileIDs] = findgroups(T.ProfileID);
TSeasonProfiles = table(profileIDs(integratedProfiles.ProfileGroup), ...
    integratedProfiles.Month, integratedProfiles.MeanChl, ...
    'VariableNames', {'ProfileID', 'Month', 'Upper30mChl'});
TSeasonProfiles.SeaIce = nan(height(TSeasonProfiles), 1);
TSeasonProfiles.Month_doy = nan(height(TSeasonProfiles), 1);
for iProfile = 1:height(TSeasonProfiles)
    id = TSeasonProfiles.ProfileID(iProfile);
    profileTimes = unique(T.(monthVar)(T.ProfileID == id));
    if numel(profileTimes) ~= 1
        error('Fram:NonUniqueProfileTime', ...
            'ProfileID %g has multiple Month_doy values; expected exactly one.', id);
    end
    TSeasonProfiles.Month_doy(iProfile) = profileTimes;
    iceValues = T.(iceVar)(T.ProfileID == id);
    iceValues = unique(iceValues(isfinite(iceValues)));
    if numel(iceValues) > 1
        error('Fram:NonUniqueProfileIce', ...
            'ProfileID %g has %d distinct sea_ice_now values; expected exactly one.', ...
            id, numel(iceValues));
    end
    if ~isempty(iceValues)
        TSeasonProfiles.SeaIce(iProfile) = iceValues;
    end
end
if any(~isfinite(TSeasonProfiles.SeaIce))
    warning('Fram:MissingProfileIce', ...
        'Sea ice is missing for %d accepted profiles; ice summaries omit these profiles.', ...
        nnz(~isfinite(TSeasonProfiles.SeaIce)));
end

%% Bin coordinates remain the mean coordinate of samples in the layer
[TSeasonChlLon, TSeasonIceLon, lonBinEdges, lonSampleRows] = local_binStatistics( ...
    tStPredictors, layerRows, TSeasonProfiles, 'Lon', lonBinWidth, profilePercentiles, ...
    seasonStart, seasonEnd, monthBinStep, minProfilesPerBin);
[TSeasonChlTemp, TSeasonIceTemp, tempBinEdges, tempSampleRows, tempBubbleCoordinates] = local_binStatistics( ...
    tStPredictors, layerRows, TSeasonProfiles, 'temp', tempBinWidth, profilePercentiles, ...
    seasonStart, seasonEnd, monthBinStep, minProfilesPerBin);

%% Four panels: chlorophyll spans two rows, sea ice occupies one row
[fg, MaxFigure] = addFindFigure('Chl_LongTempIce_Season_Upper30m', MaxFigure, ...
    [10, 20, 800, 470], 'normal', 'on');
tl = tiledlayout(fg, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
chlColor = [0.05 0.40 0.25]; % dark green, distinct from the brighter Type 5 green
iceColor = [0 0.45 0.74];

ax1 = nexttile(tl, 1, [2 1]);
local_plotPanel(ax1, TSeasonChlLon, lonBinEdges, 'Lon', ...
    'SeasonalChl', '', 'Chlorophyll-a, upper 30 m (\mug L^{-1})', chlColor);
if plotCommunityBubbles
    hold(ax1, 'on');
    lonBubbleHandles = local_bubbleChlByCommunityType(ax1, tStPredictors, tDMC, ...
        lonSampleRows, 'Lon', chlVar, lonBinEdges);
    hold(ax1, 'off');
end

ax2 = nexttile(tl, 2, [2 1]);
local_plotPanel(ax2, TSeasonChlTemp, tempBinEdges, 'temp', ...
    'SeasonalChl', '', '', chlColor);
if plotCommunityBubbles
    hold(ax2, 'on');
    typeLegendHandles = local_bubbleChlByCommunityType(ax2, tStPredictors, tDMC, ...
        tempSampleRows, 'temp', chlVar, tempBinEdges, tempBubbleCoordinates);
    hold(ax2, 'off');
end

ax3 = nexttile(tl, 5);
local_plotPanel(ax3, TSeasonIceLon, lonBinEdges, 'Lon', ...
    'SeasonalIce', 'Longitude', 'Sea ice concentration (%)', iceColor);

ax4 = nexttile(tl, 6);
local_plotPanel(ax4, TSeasonIceTemp, tempBinEdges, 'temp', ...
    'SeasonalIce', 'Average water temperature (°C)', '', iceColor);
linkaxes([ax1 ax3], 'x');
linkaxes([ax2 ax4], 'x');
ylim(ax1, [0 4]);
ylim(ax2, [0 4]);

if plotCommunityBubbles
    maxBubbleSamples = 50;
    for hBubble = [lonBubbleHandles; typeLegendHandles]'
        maxBubbleSamples = max(maxBubbleSamples, max(hBubble.SizeData));
    end
    bubbleSampleLimits = [2 maxBubbleSamples];
    bubbleDiameterLimits = sqrt(bubbleAreaPerSample * bubbleSampleLimits);
    bubblesize(ax1, bubbleDiameterLimits);
    bubblesize(ax2, bubbleDiameterLimits);
    bubblelim(ax1, bubbleSampleLimits);
    bubblelim(ax2, bubbleSampleLimits);
    if ~isempty(typeLegendHandles)
        lg = legend(ax2, typeLegendHandles, 'Location', 'northeast', 'AutoUpdate', 'off');
        lg.Title.String = 'Community type';
        lg.Layout.Tile = 'east';
    end
    if ~isempty(typeLegendHandles)
        bubbleLegendValues = [bubbleSampleLimits(1), mean(bubbleSampleLimits), bubbleSampleLimits(2)];
        lgSize = bubblelegend(ax2, 'Samples per bubble', ...
            'Location', 'northeast', 'NumBubbles', 3, 'BubbleSizeOrder', 'ascending', ...
            'LimitLabels', cellstr(compose('%g', bubbleLegendValues)));
        lgSize.Layout.Tile = 'east';
    end
end

%drawnow
AddLetters2Plots(fg, {'a', 'b', 'c', 'd'},  "Location",     "NorthEast", "FitLocation", true);


disp(TSeasonChlLon);
disp(TSeasonIceLon);
disp(TSeasonChlTemp);
disp(TSeasonIceTemp);

function [chlStats, iceStats, binEdges, sampleRows, bubbleCoordinates] = local_binStatistics( ...
    samples, layerRows, profiles, binVar, binWidth, profilePercentiles, ...
    seasonStart, seasonEnd, monthBinStep, minProfilesPerBin)
rows = layerRows & isfinite(samples.(binVar));
if ~any(rows)
    error('No finite %s observations in the seasonal 0-30 m layer.', binVar);
end
[group, ids] = findgroups(samples.ProfileID(rows));
coordinates = splitapply(@mean, samples.(binVar)(rows), group);
binMin = floor(min(coordinates) / binWidth) * binWidth;
binMax = ceil(max(coordinates) / binWidth) * binWidth;
binEdges = binMin:binWidth:binMax;
if numel(binEdges) < 2
    binEdges = [binMin, binMin + binWidth];
end
centers = binEdges(1:end-1)' + diff(binEdges)' / 2;
[hasCoordinate, loc] = ismember(profiles.ProfileID, ids);
profileCoordinates = nan(height(profiles), 1);
profileCoordinates(hasCoordinate) = coordinates(loc(hasCoordinate));
bubbleCoordinates = samples.(binVar);
if strcmp(binVar, 'temp')
    % Keep every sample of an accepted profile in that profile's temperature bin.
    [hasProfile, profileLoc] = ismember(samples.ProfileID, profiles.ProfileID);
    bubbleCoordinates = nan(height(samples), 1);
    bubbleCoordinates(hasProfile) = profileCoordinates(profileLoc(hasProfile));
    rows = layerRows & isfinite(bubbleCoordinates);
end
profileBin = zeros(height(profiles), 1);
[~, ~, profileBin(hasCoordinate)] = histcounts(profileCoordinates(hasCoordinate), binEdges);
[retainedProfiles, retainedSamples, keepBins, ~, binSampleCounts] = chl_profileBinMask( ...
    profileCoordinates, bubbleCoordinates(rows), binEdges, minProfilesPerBin);
sampleRows = false(height(samples), 1);
sampleRows(rows) = retainedSamples;
profiles = profiles(retainedProfiles, :);
profileBin = profileBin(retainedProfiles);
chlStats = local_summarize(profiles, profileBin, centers, binVar, ...
    'Upper30mChl', 'SeasonalChl', profilePercentiles, seasonStart, seasonEnd, monthBinStep);
iceStats = local_summarize(profiles, profileBin, centers, binVar, ...
    'SeaIce', 'SeasonalIce', profilePercentiles, seasonStart, seasonEnd, monthBinStep);
chlStats.NSamples = binSampleCounts;
iceStats.NSamples = binSampleCounts;
chlStats = chlStats(keepBins, :);
iceStats = iceStats(keepBins, :);
end

function summary = local_summarize(profiles, binID, centers, binVar, valueVar, seasonalVar, ...
    profilePercentiles, seasonStart, seasonEnd, monthBinStep)
nBins = numel(centers);
values = nan(nBins, 3);
counts = zeros(nBins, 3);
for b = 1:nBins
    rows = binID == b & isfinite(profiles.(valueVar));
    if ~any(rows)
        continue
    end
    v = profiles.(valueVar)(rows);
    [seasonalMean, nTimeBins, nMonths] = chl_seasonTimeBinMean( ...
        v, profiles.Month_doy(rows), seasonStart, seasonEnd, monthBinStep);
    values(b, :) = [seasonalMean, prctile(v, profilePercentiles)];
    counts(b, :) = [numel(v), nMonths, nTimeBins];
end
summary = table(centers, values(:, 1), values(:, 2), values(:, 3), ...
    counts(:, 1), counts(:, 2), counts(:, 3), 'VariableNames', ...
    {binVar, seasonalVar, 'ProfileLower', 'ProfileUpper', 'NProfiles', 'NMonths', 'NTimeBins'});
summary.Properties.UserData.ProfilePercentiles = profilePercentiles;
summary.Properties.UserData.MonthBinStep = monthBinStep;
end

function local_plotPanel(ax, stats, edges, binVar, valueVar, xLabel, yLabel, color)
hold(ax, 'on');
[~, ~, plotBins] = histcounts(stats.(binVar), edges);
for b = 1:height(stats)
    if ~isfinite(stats.ProfileLower(b)) || ~isfinite(stats.ProfileUpper(b))
        continue
    end
    j = plotBins(b);
    patch(ax, [edges(j), edges(j+1), edges(j+1), edges(j)], ...
        [stats.ProfileLower(b), stats.ProfileLower(b), stats.ProfileUpper(b), stats.ProfileUpper(b)], ...
        color, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
% Preserve line gaps at excluded bins and the original patch widths.
plotCenters = edges(1:end-1)' + diff(edges)' / 2;
plotValues = nan(size(plotCenters));
plotValues(plotBins) = stats.(valueVar);
plot(ax, plotCenters, plotValues, 's-', ...
    'LineWidth', 1.8, 'Color', color, 'MarkerFaceColor', color, ...
    'DisplayName', 'Mean of nonempty time-bin profile means', 'HandleVisibility', 'off');
hold(ax, 'off');
xlabel(ax, xLabel);
ylabel(ax, yLabel);
xlim(ax, [edges(1), edges(end)]);
end

function typeHandles = local_bubbleChlByCommunityType( ...
    ax, samples, communities, rows, binVar, chlVar, edges, bubbleCoordinates)
typeHandles = gobjects(0, 1);
if nargin < 8
    bubbleCoordinates = samples.(binVar);
end
if ~ismember('ClustID', communities.Properties.VariableNames)
    warning('Column ClustID was not found. Community-type bubbles are not plotted.');
    return
end
chl = samples.(chlVar)(rows);
types = communities.ClustID(rows);
[~, ~, bins] = histcounts(bubbleCoordinates(rows), edges);
valid = isfinite(chl) & isfinite(types) & bins > 0;
if ~any(valid)
    return
end
chl = chl(valid);
types = types(valid);
bins = bins(valid);
[groups, groupBins, groupTypes] = findgroups(bins, types);
groupN = splitapply(@numel, chl, groups);
groupMean = splitapply(@mean, chl, groups);
colors = f_Clrs_fresh(max(groupTypes));
centers = edges(1:end-1) + diff(edges) / 2;
for k = unique(groupTypes)'
    idx = groupTypes == k;
    typeHandles(end + 1, 1) = bubblechart(ax, centers(groupBins(idx)), groupMean(idx), ...
        groupN(idx), colors(k, :), 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], 'UserData', groupN(idx), ...
        'DisplayName', ['Type ' num2str(k)]);
end
end
