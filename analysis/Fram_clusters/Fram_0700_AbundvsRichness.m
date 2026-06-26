%Find abundance vs richness across phyla 

proj_Initialize;

%% Load data
[Data, DI] = dt_loadData(HGP, []);
tSamplesNorm = Data.tSamplesNorm;
dSamplesNorm = Data.dSamplesNorm;
tPhylaRoles = Data.tPhylaRoles;
tTraits = Data.tTraits;

idxSamplCols = DI.idxSamplCols;
tStPredictors = Data.tStPredictors;
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
D = load(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'));
clust_model = D.clust_model;

%% Get list of grouping taxonomic level (e.g. Phylum) 
groupTaxa = unique(tSamplesNorm(:, HGP.Samples.topGroupName));
roles = unique(tPhylaRoles.RoleName);

%find richness and relative abundance by groups phyla, diversity 
[tTaxaDiversity, tTaxaDiversitySamples] = taxaImpact(tSamplesNorm,HGP.Samples.topGroupName,'diversity',  idxSamplCols);
tTaxaDiversity = tTaxaDiversity(tPhylaRoles.Phylum, :);
tTaxaDiversitySamples = tTaxaDiversitySamples(tPhylaRoles.Phylum, :);
tTaxaDiversitySamples = (tTaxaDiversitySamples);
%tTaxaDiversitySamples = log10(tTaxaDiversitySamples./tTaxaDiversity.mean_diversity);

%Transpose data 
tTaxaDiversitySamplesTrans = rows2vars(tTaxaDiversitySamples);
tTaxaDiversitySamplesTrans.Row = tTaxaDiversitySamplesTrans.OriginalVariableNames;
tTaxaDiversitySamplesTrans.OriginalVariableNames = [];

[tTaxaAbundance, tTaxaAbundanceSamples]  = taxaImpact(tSamplesNorm,HGP.Samples.topGroupName,'abundance',  idxSamplCols);
tTaxaAbundance = tTaxaAbundance(tPhylaRoles.Phylum, :);
tTaxaAbundanceSamples = tTaxaAbundanceSamples(tPhylaRoles.Phylum, :);
tTaxaAbundanceSamples = (tTaxaAbundanceSamples).^0.25;
%tTaxaAbundanceSamples = log10(tTaxaAbundanceSamples./tTaxaAbundance.mean_diversity);

%%
[fg, MaxFigure] = addFindFigure('adund_vs_richness_taxa', MaxFigure, [-1, -1,  1218 853], 'normal', 'on');
role_colors = f_Clrs_fresh(height(tTaxaAbundance));
clearvars 'leg_handles';
for iR = 1:height(tTaxaAbundance)
    nexttile
    taxaName = tTaxaAbundance.Row(iR);
    disp(['Phylum ' taxaName])
    roleRichness = tTaxaDiversitySamples{taxaName, :};
    roleAbund = tTaxaAbundanceSamples{taxaName, :};
    s = scatter(roleRichness(:), roleAbund(:), 50, role_colors(iR, :), "filled");
    s.MarkerFaceAlpha = 0.2;
    hold on
    leg_handles(iR) = s;
    disp([taxaName])
    indValid = isfinite(roleRichness(:)) & isfinite(roleAbund(:))
    roleRichness = roleRichness(indValid);
    roleAbund = roleAbund(indValid);
    lm = fitlm(roleRichness, roleAbund)
    richnRange = linspace(min(roleRichness), max(roleRichness), 20);
    if lm.Coefficients.pValue(2)<0.01
        plot(richnRange, lm.predict(richnRange'), "Color", 'k', "LineWidth", 3)
    else
        plot(richnRange, lm.predict(richnRange'), "Color", 'k', "LineStyle","--")
    end
    title(taxaName)
    xlabel('Richness ')
    ylabel('Abundance')
end
xlim()
%legend(leg_handles, roles', 'Location','best')

%ylim([-1, 1.2])

%% plot cross correlations in richness between phyla 
%take tTaxaDiversitySamples

%find cross correlations 
[fg, MaxFigure] = addFindFigure('richnessCrossCorr', MaxFigure, [-1, -1,  1362 1362], 'normal', 'on');
R = corrplot(tTaxaDiversitySamplesTrans)
 

%%
[fg, MaxFigure] = addFindFigure('adund_vs_richness_roles', MaxFigure, [-1, -1,  1362 281], 'normal', 'on');
role_colors = f_Clrs_fresh(length(roles));
clearvars 'leg_handles';
for iR = 1:length(roles)
    nexttile
    role = roles(iR);
    ind = strcmp(tPhylaRoles.RoleName, role);
    disp(['Role ' role])
    disp(tPhylaRoles.Phylum(ind)')
    disp([tTaxaAbundance.Row(ind)']);
    disp([tTaxaDiversity.Row(ind)']);
    roleRichness = tTaxaDiversitySamples{ind, :};
    roleAbund = tTaxaAbundanceSamples{ind, :};
    s = scatter(roleRichness(:), roleAbund(:), 50, role_colors(iR, :), "filled");
    s.MarkerFaceAlpha = 0.5;
    hold on
    leg_handles(iR) = s;
    disp(['Role ' role])
    lm = fitlm(roleRichness(:), roleAbund(:))
    richnRange = -0.6:0.1:0.6;
    if lm.Coefficients.pValue(2)<0.01
        plot(richnRange, lm.predict(richnRange'), "Color", role_colors(iR, :), "LineWidth", 3)
    else
        plot(richnRange, lm.predict(richnRange'), "Color", role_colors(iR, :), "LineStyle","--")
    end
    title(role)
    xlabel('Richness ')
    ylabel('Abundance ')
end
xlim()
%legend(leg_handles, roles', 'Location','best')

ylim([-1, 1.2])
%%


impact_type = 'esn'; 
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, impact_type, idxSamplCols);
tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
HGP.Samples.topGroupNameNMostDiverse = min(HGP.Samples.topGroupNameNMostDiverse, height(tTaxaDiversity))
tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);


tTaxaDiversity = tTaxaDiversity(tPhylaRoles.Phylum, :);

[fg, MaxFigure] = addFindFigure('CommCompos divers', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_div = [];
clustStats = {};
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);

    tClustTaxaDiversity =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, impact_type, samples_idx);
    tClustTaxaDiversity = tClustTaxaDiversity(tTaxaDiversity.Row, :);
    tClustTaxaDiversity.(HGP.Samples.topGroupName) = [];
    low_y = -1; high_y = 1;
    [~, clustStats{k}] = violLogWithNames(tClustTaxaDiversity, tTaxaDiversity, clrs(k,:), [low_y,  high_y], 'YLabel','Diversity deviation index');
    title(['Type ' NS(k), ', n=' NS(length(samples_idx))]);    
end
AddLetters2Plots(fg,  'HShift', -0.02, 'VShift', 0.0, 'Location', 'NorthEast');
%D = (pdist(commmun_deviations,'correlation'));
% [fg, MaxFigure] = addFindFigure('Comm-similarities distributions', MaxFigure, [10, 20, 400 400], 'normal', 'on');
% histogram(D(:), 10)
Stat.ESNRelClust = clustStats;
Stat.ESNtaxa = tTaxaDiversity;


%% Get most rich taxa 
impact_type = 'diversity'; 
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, impact_type, idxSamplCols);
tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
HGP.Samples.topGroupNameNMostDiverse = min(HGP.Samples.topGroupNameNMostDiverse, height(tTaxaDiversity))
tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);


tTaxaDiversity = tTaxaDiversity(tPhylaRoles.Phylum, :);

[fg, MaxFigure] = addFindFigure('CommCompos richn', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_div = [];
clustStats = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);

    tClustTaxaDiversity =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, impact_type, samples_idx);
    tClustTaxaDiversity = tClustTaxaDiversity(tTaxaDiversity.Row, :);
    tClustTaxaDiversity.(HGP.Samples.topGroupName) = [];
    low_y = -1; high_y = 1;
    [~, clustStats{k}] = violLogWithNames(tClustTaxaDiversity, tTaxaDiversity, clrs(k,:), [low_y,  high_y], 'YLabel','Richness deviation index');
    title(['Type ' NS(k), ', n=' NS(length(samples_idx))]);    
end
AddLetters2Plots(fg,  'HShift', -0.02, 'VShift', 0.0, 'Location', 'NorthEast');
Stat.RichnRelClust = clustStats;
Stat.RichnTaxa = tTaxaDiversity;

%% Get most abundant phyla 
tTaxaAbundance =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'abundance', idxSamplCols);
tTaxaAbundance = tTaxaAbundance(tTaxaDiversity.Row, :);
%tTaxaAbundance = sortrows(tTaxaAbundance, "mean_diversity", "descend");
%tTaxaAbundance = tTaxaAbundance(1:HGP.Samples.topGroupNameNMostDiverse, :);

[fg, MaxFigure] = addFindFigure('CommCompos abund', MaxFigure, [-1, -1, 1.7087e+03 682], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
clustStats = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpactSamples(tSamplesNorm(:, :),HGP.Samples.topGroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    tClustTaxaAbundance.(HGP.Samples.topGroupName) = [];

    %low_y = 0.1; high_y = 10;
    low_y = -1.5; high_y = -low_y;
    [~, clustStats{k}] = violLogWithNames(tClustTaxaAbundance, tTaxaAbundance, clrs(k,:), [low_y,  high_y], 'YLabel','Abundance deviation index');
    title(['Type ' NS(k), ', n=' NS(length(samples_idx))]);    
end
AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');

Stat.AbundRelClust = clustStats;
Stat.AbundTaxa = tTaxaAbundance;

%% plot heat map of relative median values
[fg, MaxFigure] = addFindFigure('dev_heatmaps', MaxFigure, [-1, -1, 1.7087e+03 682], 'normal', 'on');
clf
cellStatistics = {Stat.AbundRelClust, Stat.RichnRelClust};  %get statistics with relative values
clims  = [-1, 1; -0.5, 0.5];
HeatMaps = {};
for iS = 1:length(cellStatistics)
    cellStats = cellStatistics{iS};
    HeatMap = [];
    for k = 1:clust_model.K
        HeatMap = [HeatMap; cellStats{k}.log_norm_val_median(:)'];
    end
    HeatMaps{iS} = HeatMap;

    nexttile
    hm = heatmap(cellStats{1}.names, 1:clust_model.K, HeatMap);
    %hm.MissingDataColor = [1, 1, 1]*0.9;
    %hm.CellLabelFormat = '%0.4g';
    colormap(redblue); caxis(clims(iS, :));
    cellStatistics{iS} = cellStats;
end
Stat.AbundHeatMap = HeatMaps{1};
Stat.RichnHeatMap = HeatMaps{2};
%% plot abundance vs richness 
[fg, MaxFigure] = addFindFigure('adund_vs_richness', MaxFigure, [-1, -1,  1362 281], 'normal', 'on');
roles = unique(tPhylaRoles.RoleName);
role_colors = f_Clrs_fresh(length(roles));
clearvars 'leg_handles';
for iR = 1:length(roles)
    nexttile
    role = roles(iR);  
    ind = strcmp(tPhylaRoles.RoleName, role);
    disp(['Role ' role])
    disp(tPhylaRoles.Phylum(ind)')
    disp([Stat.AbundRelClust{1}.names(ind)']);
    disp([Stat.RichnRelClust{1}.names(ind)']);
    roleRichness = Stat.RichnHeatMap(:, ind);
    roleAbund = Stat.AbundHeatMap(:, ind);
    s = scatter(roleRichness(:), roleAbund(:), 50, role_colors(iR, :), "filled");
    s.MarkerFaceAlpha = 0.5;
    hold on
    leg_handles(iR) = s;
    disp(['Role ' role])
    lm = fitlm(roleRichness(:), roleAbund(:))
    richnRange = -0.6:0.1:0.6;
    if lm.Coefficients.pValue(2)<0.01
        plot(richnRange, lm.predict(richnRange'), "Color", role_colors(iR, :), "LineWidth", 3)
    else
        plot(richnRange, lm.predict(richnRange'), "Color", role_colors(iR, :), "LineStyle","--")
    end
    title(role)
    xlabel('Richness deviation')
ylabel('Abundance deviation')
end
xlim()
%legend(leg_handles, roles', 'Location','best')

ylim([-1, 1.2])

  
%% Get most abundant species 

GroupName = 'Species';
%GroupName = 'Genus';
%GroupName = 'Family';
%GroupName = 'Order';
%GroupName = 'Class';

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
    tTaxaAbundance = tTaxaAbundance(1:20, :);
    [iloca, ilocb] = ismember(tTaxaAbundance.(GroupName), tSamplesNormMerged.(GroupName));
    T = [tSamplesNormMerged(ilocb, 1:7), tTaxaAbundance(:, 'mean_diversity')];
    T = renamevars(T, 'mean_diversity', 'mean_abundance');
    writenewtable(T, fullfile(HGP.Data.PathMetaData, "mainSpecies.xls"));
end 


tTaxaAbundance = tTaxaAbundance(1:20, :);
tTaxaAbundance = sortrows(tTaxaAbundance, GroupName, "ascend");


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


   %low_y = 0.1; high_y = 10;
    low_y = -2; high_y = -low_y;
    violLogWithNames(tClustTaxaAbundance, tTaxaAbundance, clrs(k,:), [low_y,  high_y], 'YLabel','Abundance deviation index')
    title(['Type ' NS(k), ', n=' NS(length(samples_idx))]);    
end
AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');





%% plot community type vs chlorophyll 
[fg, MaxFigure] = addFindFigure('comm_type_chl', MaxFigure, [10, 20, 1.1003e+03 310.6667], 'normal', 'on');
nexttile
%boxplot(tStPredictors.chlorophyll_a_0, tDMC.ClustID)

vs = violinplot(tStPredictors.chlorophyll_a_0, tDMC.ClustID, 'ViolinColor', clrs);
%yscale('log')
ylabel('Chl-a')
xlabel('Community type')

%% plot community type, chlorophyll, vs depth
[fg, MaxFigure] = addFindFigure('comm_type_chl_depth', MaxFigure, [10, 20, 1.1003e+03 310.6667], 'normal', 'on');
for iC = 1:clust_model.K
    nexttile
    %boxplot(tStPredictors.chlorophyll_a_0, tDMC.ClustID)
    ind =  tDMC.ClustID == iC;
    h = sum(ind);
    colormap(jet)
    s = scatter(tStPredictors.chlorophyll_a_0(ind), tStPredictors.Depth(ind)+2*(rand(h, 1)-0.5), 20, tStPredictors.doy(ind), "filled");
    %clim([-2, 8])
    clim([150, 280])
    ylim([0, 100])
    xlim([0, 8])
    xscale('log')
    ax = gca;
    ax.YDir ="reverse"
end
%yscale('log')
%ylabel('Depth')
%xlabel('Chl')


%% mixed model for chl 
%make and clean input table
%% --- Fit mixed model (run once) ---
T = tStPredictors(:, {'Depth','chlorophyll_a_0'});
T.ClustID = (tDMC.ClustID);
T = rmmissing(T, 'DataVariables', {'Depth','chlorophyll_a_0','ClustID'});
[chTrans, lambda] = boxcox(T.chlorophyll_a_0(~isnan(T.chlorophyll_a_0)))
%lambda = 0.2. Let's use 1/4
lambda = 0.25;
T.transfChl = T.chlorophyll_a_0.^lambda;

% Mixed model: random intercept only (based on your compare results)
lme1 = fitlme(T, 'transfChl ~ Depth + (1|ClustID)');

% Fixed-effect-only model (optional, for marginal line)
lme0 = fitlme(T, 'transfChl ~ Depth');

compare(lme0, lme1)

lme2 = fitlme(T, 'transfChl ~ Depth + (Depth|ClustID)');
compare( lme1, lme2)

% we use lme1, because lme2 does not give any improvement 
% Random intercept estimates (BLUPs) per cluster
RE = randomEffects(lme1);  % table

% Depth grid for prediction
zGrid = linspace(0, 100, 200)';

%% --- Plot ---
[fg, MaxFigure] = addFindFigure('comm_type_chl_depth', MaxFigure, [10, 20, 1.1003e+03 310.6667], 'normal', 'on');


for iC = 1:clust_model.K
    nexttile

    ind = (tDMC.ClustID == iC);
    h = sum(ind);

    % Scatter: chl vs depth (jittered), colored by doy
   
    % s = scatter(tStPredictors.Depth(ind) + 2*(rand(h,1)-0.5), tStPredictors.chlorophyll_a_0(ind), ...
    %             20, clrs(iC, :), "filled");
   % s = scatter(tStPredictors.Depth(ind) + 2*(rand(h,1)-0.5), tStPredictors.chlorophyll_a_0(ind), ...
   %              20, tStPredictors.temp(ind), "filled");
   %  clim([-2, 8])
   s = scatter(tStPredictors.Depth(ind) + 2*(rand(h,1)-0.5), tStPredictors.chlorophyll_a_0(ind), ...
                20, tStPredictors.Month_doy(ind), "filled");
   s.MarkerFaceAlpha = 0.5;
    clim([5, 10])
    xlim([0, 100])
    ylim([0.01, 10])
    set(gca,'YScale','log')
    % ax = gca; ax.YDir = "reverse";
    hold on

    % --- Mixed-model predictions ---
    % 1) Cluster-specific (conditional) prediction: includes random intercept for this ClustID
    tblC = table(zGrid, repmat(iC, numel(zGrid), 1), 'VariableNames', {'Depth','ClustID'});
    yHatC = predict(lme1, tblC, 'Conditional', true);   % predicted logChl
    chlHatC = yHatC.^(1/lambda);

    % 2) Population (marginal) prediction: fixed effects only
    tbl0 = table(zGrid, 'VariableNames', {'Depth'});
    yHat0 = predict(lme0, tbl0);                        % predicted logChl
    chlHat0 = (yHat0).^(1/lambda);

    % Plot lines (y = chl, x = depth)
    plot(zGrid, chlHat0, 'k-', 'LineWidth', 1.2);       % marginal (fixed only)
    plot(zGrid, chlHatC, 'r-', 'LineWidth', 1.6, 'Color', clrs(iC, :));       % conditional (cluster)
 colorbar
    % Add cluster random intercept estimate in title (optional)
    % RE has rows for each random effect; filter intercept for this cluster

    hold on
end

 ylabel('Chl-a')
 xlabel('Depth (m)')


