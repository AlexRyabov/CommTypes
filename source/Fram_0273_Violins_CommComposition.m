%Find composition of  communities from different clusters 
%Make violin plots 

proj_Initialize;

%% Load data
[Data, DI] = dt_loadData(HGP, []);
tSamplesNorm = Data.tSamplesNorm;
dSamplesNorm = Data.dSamplesNorm;
tPhylaRoles = Data.tPhylaRoles;
tTraits = Data.tTraits;
[roleClrs, roleShortNames] = Fram_roleColorsAndNames();

idxSamplCols = DI.idxSamplCols;
tStPredictors = Data.tStPredictors;
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
D = load(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'));
clust_model = D.clust_model;
plotLabels = local_plotLabels(clust_model.K);

%% Role color legend
[fgRoleLegend, MaxFigure] = addFindFigure('Role color legend', MaxFigure, [10, 20, 380, 180], 'normal', 'on');
clf(fgRoleLegend);
axRoleLegend = axes(fgRoleLegend);
roleOrder = {'Autotroph', 'Mixotroph', 'Zooplankton', 'Parasite', 'Recycler'};
hold(axRoleLegend, 'on');

for iRole = 1:numel(roleOrder)
    roleName = roleOrder{iRole};
    scatter(axRoleLegend, iRole, 1, 80, roleClrs(roleName), 'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'DisplayName', roleShortNames(roleName));
end

plot(axRoleLegend, 1:numel(roleOrder), ones(1, numel(roleOrder)), '-', ...
    'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
hold(axRoleLegend, 'off');
box(axRoleLegend, 'on');
grid(axRoleLegend, 'on');
xlim(axRoleLegend, [0.5 numel(roleOrder) + 0.5]);
ylim(axRoleLegend, [0.5 1.5]);
axRoleLegend.YTick = [];
axRoleLegend.XTick = 1:numel(roleOrder);
axRoleLegend.XTickLabel = roleOrder;
axRoleLegend.XTickLabelRotation = 30;
legend(axRoleLegend, 'Location', 'eastoutside');
title(axRoleLegend, 'Ecological roles');

%% Get most diverse taxa 
impact_type = 'esn'; 
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, impact_type, idxSamplCols);
tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
HGP.Samples.topGroupNameNMostDiverse = min(HGP.Samples.topGroupNameNMostDiverse, height(tTaxaDiversity))
tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);

isPresent = ismember(tPhylaRoles.Phylum, tTaxaDiversity.Row);
tPhylaRoles  = tPhylaRoles(isPresent, :);
tTaxaDiversity = tTaxaDiversity(tPhylaRoles.Phylum, :);
taxaRoleColors = local_roleColorsForNames(tTaxaDiversity.Row, tPhylaRoles, roleClrs);

[fg, MaxFigure] = addFindFigure('CommCompos divers', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_div = [];
SamplesStat = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);

    tClustTaxaDiversity =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, impact_type, samples_idx);
    tClustTaxaDiversity = tClustTaxaDiversity(tTaxaDiversity.Row, :);
    tClustTaxaDiversity.(HGP.Samples.topGroupName) = [];
    low_y = -1; high_y = 1;
    [~, SamplesStat] = violLogWithNames(tClustTaxaDiversity, tTaxaDiversity, taxaRoleColors, [low_y,  high_y], 'YLabel','Diversity deviation');
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg, plotLabels, 'HShift', -0.03, 'VShift', -0.02);

%%


%% Get most rich taxa 
impact_type = 'diversity'; 
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, impact_type, idxSamplCols);
tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
HGP.Samples.topGroupNameNMostDiverse = min(HGP.Samples.topGroupNameNMostDiverse, height(tTaxaDiversity))
tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);


tTaxaDiversity = tTaxaDiversity(tPhylaRoles.Phylum, :);
taxaRoleColors = local_roleColorsForNames(tTaxaDiversity.Row, tPhylaRoles, roleClrs);

[fg, MaxFigure] = addFindFigure('CommCompos richn', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_div = [];
SamplesStd = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);

    tClustTaxaDiversity =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, impact_type, samples_idx);
    tClustTaxaDiversity = tClustTaxaDiversity(tTaxaDiversity.Row, :);
    tClustTaxaDiversity.(HGP.Samples.topGroupName) = [];
    low_y = -1; high_y = 1;
    [~, SamplesStd] = violLogWithNames(tClustTaxaDiversity, tTaxaDiversity, taxaRoleColors, [low_y,  high_y], 'YLabel','Richness deviation');
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg, plotLabels, 'HShift', -0.03, 'VShift', -0.02);

%% Get most abundant phyla 
tTaxaAbundance =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'abundance', idxSamplCols);
tTaxaAbundance = tTaxaAbundance(tTaxaDiversity.Row, :);
taxaRoleColors = local_roleColorsForNames(tTaxaAbundance.Row, tPhylaRoles, roleClrs);

[fg, MaxFigure] = addFindFigure('CommCompos abund', MaxFigure, [-1, -1, 1127 537], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
SamplesStd = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    tClustTaxaAbundance.(HGP.Samples.topGroupName) = [];

    low_y = -1.5; high_y = -low_y;
    [~, SamplesStd] = violLogWithNames(tClustTaxaAbundance, tTaxaAbundance, taxaRoleColors, [low_y,  high_y], 'YLabel','Abundance deviation');
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg, plotLabels, 'HShift', -0.03, 'VShift', -0.02);



%% Get most abundant species 

GroupName = 'Species';

if 1 % merge Dino Group clades 
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNorm, "Species", "Dino-Group-II-", "Dino-Group-II-Clades");
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNormMerged, "Genus", "Dino-Group-II-", "Dino-Group-II-Clades");
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNormMerged, "Family", "Dino-Group-II-", "Dino-Group-II-Clades");
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNormMerged, "Species", "Dino-Group-I-", "Dino-Group-I-Clades");
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNormMerged, "Genus", "Dino-Group-I-", "Dino-Group-I-Clades");
    tSamplesNormMerged = cleandata_unifyNames(tSamplesNormMerged, "Family", "Dino-Group-I-", "Dino-Group-I-Clades");
    tSamplesNormMerged.TopName = tSamplesNormMerged.(GroupName);
    colNames = tSamplesNormMerged.Properties.VariableNames(idxSamplCols);
    aggNames = tSamplesNormMerged.Properties.VariableNames( 1:idxSamplCols(1)-1);
    tSamplesNormMergedAgg = grpstats(tSamplesNormMerged, aggNames, 'sum');
    tSamplesNormMergedAgg.Row = [];
    tSamplesNormMergedAgg.GroupCount = [];
    tSamplesNormMergedAgg.Properties.VariableNames(idxSamplCols) = colNames;
end

tTaxaAbundance =  taxaImpact(tSamplesNormMerged,GroupName, 'abundance', idxSamplCols);
tTaxaAbundance('NA', :) = [];
tTaxaAbundance = sortrows(tTaxaAbundance, "mean_diversity", "descend");
if 0 %update table with species information 
%select top 50 species to save 
    tTaxaAbundance = tTaxaAbundance(1:50, :);
    [iloca, ilocb] = ismember(tTaxaAbundance.(Genus), tSamplesNormMerged.(Genus));
    T = [tSamplesNormMerged(ilocb, 1:7), tTaxaAbundance(:, 'mean_diversity')];
    T = renamevars(T, 'mean_diversity', 'mean_abundance');
    writenewtable(T, fullfile(HGP.Data.Path2ResultsData, "mainSpecies.xls"));
end 


tTaxaAbundance = tTaxaAbundance(1:20, :);
tTaxaAbundance = sortrows(tTaxaAbundance, GroupName, "ascend");
taxaRoleColors = local_roleColorsForGroupedTaxa(tTaxaAbundance.Row, tSamplesNormMerged, GroupName, tPhylaRoles, roleClrs);


[fg, MaxFigure] = addFindFigure('CommCompos spec rel abund', MaxFigure, [10, 20,  1.7087e+03 682], 'normal', 'on');
drawnow
clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpactSamples(tSamplesNormMerged(:, :),GroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    tClustTaxaAbundance.(GroupName) = [];


    low_y = -2; high_y = -low_y;
    violLogWithNames(tClustTaxaAbundance, tTaxaAbundance, taxaRoleColors, [low_y,  high_y], 'YLabel','Abundance deviation')
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg, plotLabels, 'HShift', -0.03, 'VShift', -0.02);

%%
[fg, MaxFigure] = addFindFigure('CommCompos spec abund', MaxFigure, [10, 20,  1.7087e+03 682], 'normal', 'on');
drawnow
clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpactSamples(tSamplesNormMerged(:, :),GroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    tClustTaxaAbundance.(GroupName) = [];


    low_y = -5; high_y = 0;
    tTaxaAbundance.mean_diversity(:) = 1; 
    violLogWithNames(tClustTaxaAbundance, tTaxaAbundance, taxaRoleColors, [low_y,  high_y], 'YLabel','log_{10) Abundance', 'YScale','log')
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg, plotLabels, 'HShift', -0.03, 'VShift', -0.02);


%% plot community type vs chlorophyll 
[fg, MaxFigure] = addFindFigure('comm_type_chl', MaxFigure, [10, 20, 1100 175], 'normal', 'on');
nexttile

vs = violinplot(tStPredictors.chlorophyll_a_0, tDMC.ClustID, 'ViolinColor', clrs);
ylabel('Chl-a (\mug L^{-1})')
xlabel('Community type')


function labels = local_plotLabels(nLabels)
alphabet = 'abcdefghijklmnopqrstuvwxyz';

if nLabels <= numel(alphabet)
    labels = cellstr(alphabet(1:nLabels)');
else
    labels = cell(nLabels, 1);

    for iLabel = 1:nLabels
        labels{iLabel} = char("a" + string(iLabel));
    end
end
end

function colors = local_roleColorsForNames(taxonNames, tPhylaRoles, roleClrs)
colors = repmat(roleClrs('Unknown'), numel(taxonNames), 1);

for iTaxon = 1:numel(taxonNames)
    idxRole = string(tPhylaRoles.Phylum) == string(taxonNames(iTaxon));

    if any(idxRole)
        roleName = local_firstRoleName(tPhylaRoles, idxRole);
        colors(iTaxon, :) = local_roleColor(roleName, roleClrs);
    end
end
end

function colors = local_roleColorsForGroupedTaxa(taxonNames, sourceTable, groupName, tPhylaRoles, roleClrs)
colors = repmat(roleClrs('Unknown'), numel(taxonNames), 1);

if ~ismember(groupName, sourceTable.Properties.VariableNames) || ~ismember('Phylum', sourceTable.Properties.VariableNames)
    return
end

for iTaxon = 1:numel(taxonNames)
    idxTaxon = string(sourceTable.(groupName)) == string(taxonNames(iTaxon));

    if any(idxTaxon)
        phylumName = string(sourceTable.Phylum(find(idxTaxon, 1, 'first')));
        idxRole = string(tPhylaRoles.Phylum) == phylumName;

        if any(idxRole)
            roleName = local_firstRoleName(tPhylaRoles, idxRole);
            colors(iTaxon, :) = local_roleColor(roleName, roleClrs);
        end
    end
end
end

function roleName = local_firstRoleName(tPhylaRoles, idxRole)
if ismember('RoleName', tPhylaRoles.Properties.VariableNames)
    roleName = string(tPhylaRoles.RoleName(find(idxRole, 1, 'first')));
elseif ismember('Role', tPhylaRoles.Properties.VariableNames)
    roleName = string(tPhylaRoles.Role(find(idxRole, 1, 'first')));
else
    roleName = "Unknown";
end
end

function clr = local_roleColor(roleName, roleClrs)
roleName = char(string(roleName));

if isKey(roleClrs, roleName)
    clr = roleClrs(roleName);
else
    clr = roleClrs('Unknown');
end
end

