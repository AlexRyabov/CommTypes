%
% Compare community types by coherent taxonomic differences
%% Methods 
% We  select taxa which have 75% prevalence in the community type and their
% mean abundance is at least 0.01. 
% Panel C:
% Bubble heatmap of characteristic taxa across community types
% Color = enrichment relative to overall mean abundance
% Dot size = prevalence within community type
%
% Panel D:
% Pairwise succession contrasts
% Contrast = log10((mean abundance * prevalence in later type) /
%                  (mean abundance * prevalence in earlier type))
%
% Requirements:
% Data.tSamplesNorm
% Data.tStPredictors
% Data.tPhylaRoles
% DI.idxSamplCols
%
% tSamplesNorm:
% first columns contain taxonomy, sample columns start at DI.idxSamplCols.
% sample columns must be in the same order as rows of tStPredictors and tDMC.
%
% tDMC:
% table with ClustID for each sample, same row order as tStPredictors.

clearvars
proj_Initialize;

%% Load data

[Data, DI] = dt_loadData(HGP, []);

tSamplesNorm  = Data.tSamplesNorm;
tStPredictors = Data.tStPredictors;
tPhylaRoles   = Data.tPhylaRoles;

idxSamplCols = DI.idxSamplCols;

% Cluster assignment table
% Change this line if your tDMC is stored elsewhere
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));

if ~ismember('ClustID', tDMC.Properties.VariableNames)
    error('Column ClustID was not found in tDMC.');
end

communityTypes = 1:8;
presentTypes = unique(tDMC.ClustID(:))';

missingTypes = setdiff(communityTypes, presentTypes);
extraTypes = setdiff(presentTypes, communityTypes);

if ~isempty(missingTypes)
    error('Missing community type(s) in tDMC.ClustID: %s.', num2str(missingTypes));
end

if ~isempty(extraTypes)
    error('Unexpected community type(s) in tDMC.ClustID: %s.', num2str(extraTypes));
end

K = numel(communityTypes);

%% Parameters

% Taxonomic level for visualization
GroupName = 'Species';

% Use only active phyla from Data.tPhylaRoles
useOnlyActivePhyla = true;

% Detection threshold for prevalence
% A taxon is treated as present if abundance > detectionThr
detectionThr = 0;

% Small value to avoid log(0)
epsVal = 1e-8;

% Panel C settings
minMaxPrevalence = 0.75;
minMaxMeanAbund  = 1e-2;
bubbleColorMode = 'MeanAbund'; % 'Enrichment' or 'MeanAbund'
enrichCLim = 1;

% Panel D settings
pairs = [1 2; 2 3;2 6;
         3 4;
         4 5;
          7 8; ];

topN_pair = 20;
minPairPrevalence = 0.75;
minPairMeanAbund  = 1e-2;
nBootPair = 500;
pairCiPrct = [2.5 97.5];

%% Prepare phylum table

if useOnlyActivePhyla
    tPhylaRoles = tPhylaRoles(tPhylaRoles.Active == 1, :);
end

if ismember('Oder', tPhylaRoles.Properties.VariableNames)
    tPhylaRoles = sortrows(tPhylaRoles, 'Oder');
else
    warning('Column Oder was not found in tPhylaRoles. Phyla will keep original order.');
end

phylaOrdered = string(tPhylaRoles.Phylum);
phylaOrdered = strip(phylaOrdered);

%% Optional: simplify some taxonomic names

if true
    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Species", "Dino-Group-II-", "Dino-Group-II-Clades");
    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Genus",   "Dino-Group-II-", "Dino-Group-II-Clades");
    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Family",  "Dino-Group-II-", "Dino-Group-II-Clades");

    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Species", "Dino-Group-I-", "Dino-Group-I-Clades");
    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Genus",   "Dino-Group-I-", "Dino-Group-I-Clades");
    tSamplesNorm = cleandata_unifyNames(tSamplesNorm, "Family",  "Dino-Group-I-", "Dino-Group-I-Clades");
end

%% Aggregate taxa by Phylum + selected taxonomic level

[taxaNames, displayNames, plotNames, X, taxInfo] = local_aggregateTaxa( ...
    tSamplesNorm, idxSamplCols, GroupName, phylaOrdered);

% X is taxa x samples
nTaxa = size(X, 1);
nSamples = size(X, 2);

if height(tDMC) ~= nSamples
    error('Number of rows in tDMC does not match number of sample columns in tSamplesNorm.');
end

if height(tStPredictors) ~= nSamples
    warning('Number of rows in tStPredictors does not match number of sample columns in tSamplesNorm.');
end

%% Add role information to taxInfo

taxInfo.RoleName = strings(height(taxInfo), 1);
taxInfo.Role     = strings(height(taxInfo), 1);

for i = 1:height(taxInfo)
    idxRole = string(tPhylaRoles.Phylum) == string(taxInfo.Phylum(i));

    if any(idxRole)
        idxFirst = find(idxRole, 1, 'first');

        taxInfo.RoleName(i) = string(tPhylaRoles.RoleName{idxFirst});
        taxInfo.Role(i)     = string(tPhylaRoles.Role{idxFirst});
    else
        taxInfo.RoleName(i) = "Unknown";
        taxInfo.Role(i)     = "u";
    end
end

autotrophTaxaOverride = [
    "Bathycoccus_prasinos"
    "Chloroparvula_pacifica"
    "Micromonas_commoda_A2"
    "Micromonas_polaris"
    "Phaeocystis_pouchetii"
];

idxAutotrophOverride = ismember(string(taxInfo.TaxonName), autotrophTaxaOverride);
taxInfo.RoleName(idxAutotrophOverride) = "Autotroph";
taxInfo.Role(idxAutotrophOverride) = "a";

%% Colors

if exist('f_Clrs_fresh', 'file')
    clrs = f_Clrs_fresh(K);
else
    clrs = lines(K);
end

% Colors for functional roles
[roleClrs, roleShortNames] = Fram_roleColorsAndNames();

%% Calculate mean abundance and prevalence per community type

meanAbund = nan(nTaxa, K);
prevalence = nan(nTaxa, K);
nSamplesType = nan(K, 1);

for k = 1:K
    typeID = communityTypes(k);
    idx = tDMC.ClustID == typeID;
    nSamplesType(k) = sum(idx);

    meanAbund(:, k) = mean(X(:, idx), 2, 'omitnan');
    prevalence(:, k) = mean(X(:, idx) > detectionThr, 2, 'omitnan');
end

overallMeanAbund = mean(X, 2, 'omitnan');

% Enrichment relative to the average abundance across all samples
enrichment = log10((meanAbund + epsVal) ./ (overallMeanAbund + epsVal));

%% Select taxa for Panel C

okTaxa = max(prevalence, [], 2) >= minMaxPrevalence & ...
         max(meanAbund, [], 2) >= minMaxMeanAbund & ...
         ~ismissing(taxaNames) & taxaNames ~= "" & taxaNames ~= "NA";

nAvailable = sum(okTaxa);

if nAvailable == 0
    error('No taxa passed the Panel C filters. Reduce minMaxPrevalence or minMaxMeanAbund.');
end

selTaxa = find(okTaxa);

% Display selected taxa in biological order:
% phylum order from tPhylaRoles, then taxon name
tSel = taxInfo(selTaxa, :);
[~, srt] = sortrows(tSel, {'PhylumOrder', 'TaxonName'});
selTaxa = selTaxa(srt);

%% Panel C: bubble heatmap


[fgC, MaxFigure] = addFindFigure('Main comp', MaxFigure, [10, 20, 1063 833], 'normal', 'on');
clf(fgC);
fgC.Color = 'w';

axC = axes(fgC);
hold(axC, 'on');

% Larger left margin for long taxon names
axC.Position = [0.26 0.10 0.54 0.82];

E = enrichment(selTaxa, :);
M = meanAbund(selTaxa, :);
P = prevalence(selTaxa, :);

switch lower(string(bubbleColorMode))
    case "enrichment"
        bubbleColor = E;
        bubbleColorLabel = 'log_{10}(mean abundance in type / overall mean abundance)';
        bubbleColorMap = local_blueWhiteRed(256);
        bubbleCLim = [-enrichCLim enrichCLim];
        bubbleColorScale = 'linear';

    case "meanabund"
        bubbleCLim = [minMaxMeanAbund / 10 1];
        bubbleColor = max(M, bubbleCLim(1));
        bubbleColorLabel = 'Mean relative abundance in community type';
        bubbleColorMap = turbo(256);
        bubbleColorScale = 'log';

    otherwise
        error('Unknown bubbleColorMode: %s. Use ''Enrichment'' or ''MeanAbund''.', bubbleColorMode);
end

[yGrid, xGrid] = ndgrid(1:numel(selTaxa), 1:K);
dotSize = 15 + 260 * P(:);

scatter(axC, xGrid(:), yGrid(:), dotSize, bubbleColor(:), ...
    'filled', ...
    'MarkerEdgeColor', [0.25 0.25 0.25], ...
    'LineWidth', 0.25);

set(axC, 'YDir', 'reverse');

xlim(axC, [0.5 K + 1.25]);
ylim(axC, [0.5 numel(selTaxa) + 0.5]);

xticks(axC, 1:K);
xticklabels(axC, compose('Type %d', communityTypes));

yticks(axC, 1:numel(selTaxa));
yticklabels(axC, local_coloredYTickLabels(plotNames(selTaxa), taxInfo.RoleName(selTaxa), roleClrs));

axC.TickLabelInterpreter = 'tex';
axC.FontSize = 8;

xlabel(axC, 'Community type');
ylabel(axC, 'Ta');

local_addRoleLegend(axC, roleClrs, roleShortNames, 8);

colormap(axC, bubbleColorMap);
axC.ColorScale = bubbleColorScale;
cb = colorbar(axC);
cb.Label.String = bubbleColorLabel;
clim(axC, bubbleCLim);

% Move colorbar slightly to the right
cb.Position = [0.84 0.12 0.018 0.78];

% No grid lines
grid(axC, 'off');
box(axC, 'on');

% Add horizontal separators only between phyla
phylaC = string(taxInfo.Phylum(selTaxa));

for ii = 1:(numel(selTaxa)-1)
    if phylaC(ii) ~= phylaC(ii+1)
        yline(axC, ii + 0.5, '-', ...
            'Color', [0.25 0.25 0.25], ...
            'LineWidth', 0.9);
    end
end

% Add phylum names on the right side only
uniquePhylaC = unique(phylaC, 'stable');

for iP = 1:numel(uniquePhylaC)
    idxP = find(phylaC == uniquePhylaC(iP));
    yMid = mean(idxP);

    text(axC, K + 0.23, yMid, uniquePhylaC(iP), ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

% Add sample numbers above columns
yl = ylim(axC);


% Prevalence legend outside the main axes
axLeg = axes(fgC, 'Position', [0.88 0.12 0.09 0.22]);
hold(axLeg, 'on');

prevVals = [0.25 0.50 0.75 1.00];
sizes = 15 + 260 * prevVals;

for i = 1:numel(prevVals)
    scatter(axLeg, 1, i, sizes(i), ...
        [0.75 0.75 0.75], ...
        'filled', ...
        'MarkerEdgeColor', [0.25 0.25 0.25], ...
        'LineWidth', 0.25);

    text(axLeg, 1.35, i, num2str(prevVals(i), '%.2f'), ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8);
end

title(axLeg, 'Prevalence', 'FontSize', 9);
xlim(axLeg, [0.7 2.2]);
ylim(axLeg, [0.5 numel(prevVals) + 0.5]);
set(axLeg, 'YDir', 'reverse');
axis(axLeg, 'off');

%% Panel C table in workspace

tBubble = taxInfo(selTaxa, :);

for k = 1:K
    tBubble.(['MeanAbund_Type' num2str(k)]) = meanAbund(selTaxa, k);
    tBubble.(['Prevalence_Type' num2str(k)]) = prevalence(selTaxa, k);
    tBubble.(['Enrichment_Type' num2str(k)]) = enrichment(selTaxa, k);
end

%% Save Panel C species list

coreTaxonVarName = matlab.lang.makeValidName(GroupName);

tCoreSpecies = table( ...
    string(taxInfo.Phylum(selTaxa)), ...
    string(taxInfo.TaxonName(selTaxa)), ...
    string(taxInfo.RoleName(selTaxa)), ...
    'VariableNames', {'Phylum', coreTaxonVarName, 'Role'});

for k = 1:K
    typeID = communityTypes(k);
    tCoreSpecies.(['MeanRelAbund_Type' num2str(typeID)]) = meanAbund(selTaxa, k);
end

writetable(tCoreSpecies, fullfile(HGP.Data.Path2ResultsData, 'CommTypeCoreSpecies.csv'));

coreSpeciesData = struct();
coreSpeciesData.GroupName = GroupName;
coreSpeciesData.Species = string(taxInfo.TaxonName(selTaxa));
coreSpeciesData.PlotName = string(taxInfo.PlotName(selTaxa));
coreSpeciesData.Phylum = string(taxInfo.Phylum(selTaxa));
coreSpeciesData.Role = string(taxInfo.RoleName(selTaxa));
coreSpeciesData.Sample = string(tSamplesNorm.Properties.VariableNames(idxSamplCols))';
coreSpeciesData.Abundance = X(selTaxa, :)';
coreSpeciesData.CommunityTypes = communityTypes(:);
coreSpeciesData.ClustID = tDMC.ClustID(:);
coreSpeciesData.MeanRelAbund = meanAbund(selTaxa, :);
coreSpeciesData.Prevalence = prevalence(selTaxa, :);

save(fullfile(HGP.Data.Path2ResultsData, 'CommTypeCoreSpeciesData.mat'), 'coreSpeciesData');

%% Panel D: pairwise succession contrasts

[fgD , MaxFigure] = addFindFigure('Main change', MaxFigure, [10, 20,1500 850], 'normal', 'on');
clf(fgD);
fgD.Color = 'w';



nPairPanels = size(pairs, 1);

tlD = tiledlayout(fgD, 'flow', 'Padding', 'compact', 'TileSpacing', 'compact');

allPairTables = table();

for iPair = 1:size(pairs, 1)

    typeA = pairs(iPair, 1);
    typeB = pairs(iPair, 2);

    if ~all(ismember([typeA typeB], communityTypes))
        warning('Skipping Type %d -> Type %d because these types are not in communityTypes.', typeA, typeB);
        continue
    end

    % Later / earlier contrast
    contrast = log10(((meanAbund(:, typeB) .* prevalence(:, typeB)) + epsVal) ./ ...
                     ((meanAbund(:, typeA) .* prevalence(:, typeA)) + epsVal));

    okPair = max(prevalence(:, [typeA typeB]), [], 2) >= minPairPrevalence & ...
             max(meanAbund(:, [typeA typeB]), [], 2) >= minPairMeanAbund & ...
             isfinite(contrast) & ...
             ~ismissing(taxaNames) & taxaNames ~= "" & taxaNames ~= "NA";

    pairScore = abs(contrast);
    pairScore(~okPair) = -inf;

    nAvailablePair = sum(isfinite(pairScore));

    if nAvailablePair == 0
        warning('No taxa passed filters for Type %d -> Type %d.', typeA, typeB);
        continue
    end

    [~, ordPair] = sort(pairScore, 'descend');
    selPair = ordPair(1:min(topN_pair, nAvailablePair));

    % Sort by phylum order, then taxon name
    tSelPair = taxInfo(selPair, :);
    [~, srtPair] = sortrows(tSelPair, {'PhylumOrder', 'TaxonName'});
    selPair = selPair(srtPair);

    contrastPlot = contrast(selPair);
    y = 1:numel(selPair);

    idxA = tDMC.ClustID == typeA;
    idxB = tDMC.ClustID == typeB;
    [ciLow, ciHigh] = local_pairContrastBootstrapCI( ...
        X, selPair, idxA, idxB, detectionThr, epsVal, nBootPair, pairCiPrct);

    axD = nexttile(tlD);
    hold(axD, 'on');

    xline(axD, 0, 'k-', 'LineWidth', 0.8);

    xMax = max(abs([contrastPlot; ciLow; ciHigh]), [], 'omitnan');
    if isempty(xMax) || xMax == 0 || isnan(xMax)
        xMax = 1;
    end

    xlim(axD, [-1.35*xMax, 1.55*xMax]);
    ylim(axD, [0.5 numel(selPair) + 0.5]);
    set(axD, 'YDir', 'reverse');

    yticks(axD, y);
    yticklabels(axD, local_coloredYTickLabels(plotNames(selPair), taxInfo.RoleName(selPair), roleClrs));
    axD.TickLabelInterpreter = 'tex';
    axD.FontSize = 7.5;

    xlabel(axD, ['log_{10}((mean abundance \times prevalence) Type ' ...
        num2str(typeB) ' / Type ' num2str(typeA) ')']);

    title(axD, ['Type ' num2str(typeA) ' \rightarrow Type ' num2str(typeB)]);

    % No grid lines
    grid(axD, 'off');
    box(axD, 'on');

    % Add confidence intervals and prevalence-scaled contrast markers
    prevEnriched = nan(numel(selPair), 1);
    pointClrs = nan(numel(selPair), 3);

    for ii = 1:numel(selPair)
        if contrastPlot(ii) >= 0
            prevEnriched(ii) = prevalence(selPair(ii), typeB);
            pointClrs(ii, :) = clrs(typeB, :);
        else
            prevEnriched(ii) = prevalence(selPair(ii), typeA);
            pointClrs(ii, :) = clrs(typeA, :);
        end

        plot(axD, [ciLow(ii) ciHigh(ii)], [y(ii) y(ii)], '-', ...
            'Color', [0.30 0.30 0.30], ...
            'LineWidth', 1.1);
    end

    scatter(axD, contrastPlot, y, 30 + 120 * prevEnriched, ...
        pointClrs, 'filled', ...
        'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeColor', [0.20 0.20 0.20], ...
        'LineWidth', 0.25);

    % Add horizontal separators only between phyla
    phylaD = string(taxInfo.Phylum(selPair));

    for ii = 1:(numel(selPair)-1)
        if phylaD(ii) ~= phylaD(ii+1)
            yline(axD, ii + 0.5, '-', ...
                'Color', [0.25 0.25 0.25], ...
                'LineWidth', 0.8);
        end
    end

    % Add phylum names on the right side only
    uniquePhylaD = unique(phylaD, 'stable');
    xl = xlim(axD);

    for iP = 1:numel(uniquePhylaD)
        idxP = find(phylaD == uniquePhylaD(iP));
        yMid = mean(idxP);

        text(axD, xl(2) - 0.02 * range(xl), yMid, uniquePhylaD(iP), ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 7, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    % Collect pair table
    tPair = taxInfo(selPair, :);
    tPair.Pair = repmat(string(['Type' num2str(typeA) '_vs_Type' num2str(typeB)]), numel(selPair), 1);
    tPair.TypeA = repmat(typeA, numel(selPair), 1);
    tPair.TypeB = repmat(typeB, numel(selPair), 1);
    tPair.Contrast_LaterOverEarlier = contrast(selPair);
    tPair.Contrast_CILow = ciLow;
    tPair.Contrast_CIHigh = ciHigh;
    tPair.MeanAbund_A = meanAbund(selPair, typeA);
    tPair.MeanAbund_B = meanAbund(selPair, typeB);
    tPair.Prevalence_A = prevalence(selPair, typeA);
    tPair.Prevalence_B = prevalence(selPair, typeB);
    tPair.Score = abs(contrast(selPair));

    if isempty(allPairTables)
        allPairTables = tPair;
    else
        allPairTables = [allPairTables; tPair]; %#ok<AGROW>
    end
end

%% Panel D table in workspace

if isempty(allPairTables)
    warning('No Panel D contrasts passed the filters.');
end

disp('Finished.');


%% Local functions

function [taxaNames, displayNames, plotNames, Xagg, taxInfo] = local_aggregateTaxa(T, idxSamplCols, GroupName, phylaOrdered)

    % Extract abundance matrix
    X = table2array(T(:, idxSamplCols));
    X = double(X);

    % Extract phylum and selected taxonomic level
    phylumRaw = string(T.Phylum);
    phylumRaw = strip(phylumRaw);

    if ~ismember(GroupName, T.Properties.VariableNames)
        error('Column %s was not found in tSamplesNorm.', GroupName);
    end

    taxaRaw = string(T.(GroupName));
    taxaRaw = strip(taxaRaw);

    % Keep only phyla present in phylaOrdered
    idxKeep = ismember(phylumRaw, phylaOrdered);

    phylumRaw = phylumRaw(idxKeep);
    taxaRaw   = taxaRaw(idxKeep);
    X         = X(idxKeep, :);

    % Remove missing or empty taxa names
    idxGood = ~(ismissing(taxaRaw) | taxaRaw == "" | taxaRaw == "NA" | taxaRaw == "NaN" | taxaRaw == "<missing>");
    phylumRaw = phylumRaw(idxGood);
    taxaRaw   = taxaRaw(idxGood);
    X         = X(idxGood, :);

    % Aggregate by Phylum + TaxonName
    key = phylumRaw + "||" + taxaRaw;
    keyUnique = unique(key, 'stable');

    nOut = numel(keyUnique);
    nSamples = size(X, 2);

    Xagg = zeros(nOut, nSamples);
    phylumOut = strings(nOut, 1);
    taxaOut   = strings(nOut, 1);

    for i = 1:nOut
        idx = key == keyUnique(i);

        Xagg(i, :) = sum(X(idx, :), 1, 'omitnan');

        firstIdx = find(idx, 1, 'first');
        phylumOut(i) = phylumRaw(firstIdx);
        taxaOut(i)   = taxaRaw(firstIdx);
    end

    % Determine phylum order from tPhylaRoles
    phylumOrder = nan(nOut, 1);

    for i = 1:nOut
        phylumOrder(i) = find(phylaOrdered == phylumOut(i), 1, 'first');
    end

    displayNames = phylumOut + "-" + taxaOut;
    plotNames = taxaOut;

    taxInfo = table(phylumOrder, phylumOut, taxaOut, displayNames, plotNames, ...
        'VariableNames', {'PhylumOrder','Phylum','TaxonName','DisplayName','PlotName'});

    % Sort by phylum order and taxon name
    [taxInfo, sortIdx] = sortrows(taxInfo, {'PhylumOrder','TaxonName'});
    Xagg = Xagg(sortIdx, :);

    taxaNames    = taxInfo.TaxonName;
    displayNames = taxInfo.DisplayName;
    plotNames    = taxInfo.PlotName;
end

function clr = local_roleColor(roleName, roleClrs)

    roleName = char(string(roleName));

    if isKey(roleClrs, roleName)
        clr = roleClrs(roleName);
    else
        clr = roleClrs('Unknown');
    end
end

function tickLabels = local_coloredYTickLabels(labels, roleNames, roleClrs)

    labels = string(labels);
    roleNames = string(roleNames);
    tickLabels = cell(numel(labels), 1);

    for i = 1:numel(labels)
        clr = local_roleColor(roleNames(i), roleClrs);
        label = local_escapeTex(labels(i));

        tickLabels{i} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%s', ...
            clr(1), clr(2), clr(3), label);
    end
end

function label = local_escapeTex(label)

    label = char(string(label));
    label = strrep(label, '\', '\\');
    label = strrep(label, '_', '\_');
    label = strrep(label, '{', '\{');
    label = strrep(label, '}', '\}');
    label = strrep(label, '^', '\^{}');
end

function [ciLow, ciHigh] = local_pairContrastBootstrapCI(X, selTaxa, idxA, idxB, detectionThr, epsVal, nBoot, ciPrct)

    XA = X(selTaxa, idxA);
    XB = X(selTaxa, idxB);

    nTaxa = numel(selTaxa);
    nA = size(XA, 2);
    nB = size(XB, 2);

    if nA == 0 || nB == 0
        ciLow = nan(nTaxa, 1);
        ciHigh = nan(nTaxa, 1);
        return
    end

    bootContrast = nan(nTaxa, nBoot);

    for iBoot = 1:nBoot
        iA = randi(nA, 1, nA);
        iB = randi(nB, 1, nB);

        XAboot = XA(:, iA);
        XBboot = XB(:, iB);

        meanA = mean(XAboot, 2, 'omitnan');
        meanB = mean(XBboot, 2, 'omitnan');
        prevA = mean(XAboot > detectionThr, 2, 'omitnan');
        prevB = mean(XBboot > detectionThr, 2, 'omitnan');

        bootContrast(:, iBoot) = log10(((meanB .* prevB) + epsVal) ./ ...
                                       ((meanA .* prevA) + epsVal));
    end

    ciLow = prctile(bootContrast', ciPrct(1))';
    ciHigh = prctile(bootContrast', ciPrct(2))';
end

function local_addRoleLegend(ax, roleClrs, roleShortNames, fontSize)

    roleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};

    x0 = -0.13;
    y0 = 1.055;
    dx = 0.12;

    for i = 1:numel(roleOrder)
        roleName = roleOrder{i};
        clr = roleClrs(roleName);
        label = roleShortNames(roleName);
        x = x0 + (i - 1) * dx;

        text(ax, x, y0, label, ...
            'Units', 'normalized', ...
            'Color', clr, ...
            'FontSize', fontSize, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
    end
end

function cmap = local_blueWhiteRed(n)

    if nargin < 1
        n = 256;
    end

    n1 = floor(n / 2);
    n2 = n - n1;

    blue  = [0.2 0.35 0.8];
    white = [1 1 1];
    red   = [0.8 0.2 0.15];

    cmap1 = [linspace(blue(1), white(1), n1)', ...
             linspace(blue(2), white(2), n1)', ...
             linspace(blue(3), white(3), n1)'];

    cmap2 = [linspace(white(1), red(1), n2)', ...
             linspace(white(2), red(2), n2)', ...
             linspace(white(3), red(3), n2)'];

    cmap = [cmap1; cmap2];
end
