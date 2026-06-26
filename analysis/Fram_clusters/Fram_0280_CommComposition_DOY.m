%dynamics of normalized phyla diversity with doy. 
%correlations between phyla diversities

proj_Initialize;

%% Load data
[Data, DI] = dt_loadData(HGP, []);
tSamplesNorm = Data.tSamplesNorm;
dSamplesNorm = Data.dSamplesNorm;
tTraits = Data.tTraits;

idxSamplCols = DI.idxSamplCols;
tStPredictors = Data.tStPredictors;
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
D = load(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'));
clust_model = D.clust_model;


%% plot richness as a function of the doy
[fg, MaxFigure] = addFindFigure('Richness on doy', MaxFigure, [10, 20, 704.6667 581.3333], 'normal', 'on');

scatter(tStPredictors.doy, tStPredictors.divRichness);
lm = fitlm(tStPredictors.doy, tStPredictors.divRichness);
hold on
x = 140:280;
hold on
plot(x, lm.predict(x'))
hold off
xlabel('DOY')
ylabel('Richness')

%% most diverse phyla 
[fg, MaxFigure] = addFindFigure('Richness phyl on doy', MaxFigure, [10, 20, 704.6667 581.3333], 'normal', 'on');
max_phyla = 5;
clrs = f_Clrs_fresh(max_phyla);
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'diversity', idxSamplCols);
for i = 1:max_phyla
    phylum = tTaxaDiversity.Phylum(i);
    ind = strcmp(tSamplesNorm.Phylum, phylum);
    phyl_asv = tSamplesNorm(ind, idxSamplCols);
    Richness = sum(phyl_asv{:, :}>0, 1);

    s = scatter(tStPredictors.doy, Richness, 'filled', 'MarkerFaceColor',clrs(i, :));
    s.Annotation.LegendInformation.IconDisplayStyle = 'off';
    s.MarkerFaceAlpha = 0.1;
    lm = fitlm(tStPredictors.doy, Richness, 'quadratic');
    hold on
    x = 140:280;
    p = plot(x, lm.predict(x'), "Color", clrs(i, :), 'LineWidth', 4);
    xlabel('DOY')
    ylabel('Richness')

end
hold off
legend( tTaxaDiversity.Phylum(1:max_phyla))

%% 
%% most diverse phyla 
[fg, MaxFigure] = addFindFigure('norm Richness phyl on doy', MaxFigure, [10, 20, 704.6667 581.3333], 'normal', 'on');
max_phyla = 10;
clrs = f_Clrs_fresh(max_phyla);
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'diversity', idxSamplCols);
tRichn = table;
for i = 1:max_phyla
    phylum = tTaxaDiversity.Phylum(i);
    ind = strcmp(tSamplesNorm.Phylum, phylum);
    phyl_asv = tSamplesNorm(ind, idxSamplCols);
    Richness = sum(phyl_asv{:, :}>0, 1);
    Richness = Richness/mean(Richness);
    tRichn.(phylum) = Richness';
    s = scatter(tStPredictors.doy, Richness, 'filled', 'MarkerFaceColor',clrs(i, :));
    s.Annotation.LegendInformation.IconDisplayStyle = 'off';
    s.MarkerFaceAlpha = 0.1;
    lm = fitlm(tStPredictors.doy, Richness, 'quadratic');
    hold on
    x = 140:280;
    p = plot(x, lm.predict(x'), "Color", clrs(i, :), 'LineWidth', 4);
    xlabel('DOY')
    ylabel('Richness')

end
f_Set_scales('logy')
maxy = 2;
ylim([1/maxy , maxy ])
%scatter(tStPredictors.doy, geomean(tRichn{:, :}, 2))
hold off
legend( tTaxaDiversity.Phylum(1:max_phyla), 'Location','southeast')

%% Ocrophytes vs chlorophytes etc
[fg, MaxFigure] = addFindFigure('correl of phyl', MaxFigure, [10, 20, 1.5233e+03 581.3333], 'normal', 'on');
for i=1:size(tRichn, 2)-1
    for j=i+1:size(tRichn, 2)
        nexttile
        s = scatter(tRichn{:, i}, tRichn{:, j}, 'filled');
        s.MarkerFaceAlpha = 0.2;
        xlabel(tRichn.Properties.VariableNames(i));
        ylabel(tRichn.Properties.VariableNames(j));
    end
end

%% Ocrophytes vs chlorophytes etc
[fg, MaxFigure] = addFindFigure('heat map', MaxFigure, [10, 20, 729.3333 581], 'normal', 'on');
[c, rho]  = corr(tRichn{:, :});
c(rho > 0.01) = NaN;

heatmap(tRichn.Properties.VariableNames, tRichn.Properties.VariableNames, c)


%%
[fg, MaxFigure] = addFindFigure('mean std ', MaxFigure, [10, 20, 400, 400], 'normal', 'on');
indPhylas = 1:max_phyla;
tRichn.mean(:) = mean(tRichn{:, indPhylas}, 2);
tRichn.std(:) = std(tRichn{:, indPhylas}, 0, 2);
scatter(tRichn.mean, tRichn.std, 20,  tStPredictors.Depth,  'filled');
colorbar
xlabel('mean');
ylabel('std')

%%
[fg, MaxFigure] = addFindFigure('Richness phyl on doy only regr', MaxFigure, [10, 20, 704.6667 581.3333], 'normal', 'on');
max_phyla = 5;
clrs = f_Clrs_fresh(max_phyla);
tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'diversity', idxSamplCols);
for i = 1:max_phyla
    phylum = tTaxaDiversity.Phylum(i);
    ind = strcmp(tSamplesNorm.Phylum, phylum);
    phyl_asv = tSamplesNorm(ind, idxSamplCols);
    Richness = sum(phyl_asv{:, :}>0, 1);

    s = scatter(tStPredictors.doy, Richness, 'filled', 'MarkerFaceColor',clrs(i, :));
    s.Annotation.LegendInformation.IconDisplayStyle = 'off';
    lm = fitlm(tStPredictors.doy, Richness, 'quadratic');
    hold on
    x = 140:280;
    p = plot(x, lm.predict(x'), "Color", clrs(i, :), 'LineWidth', 4);
    xlabel('DOY')
    ylabel('Richness')

end
hold off
legend( tTaxaDiversity.Phylum(1:max_phyla))


%% Get most diverse taxa 




tTaxaDiversity =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'diversity', idxSamplCols);
tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);

[fg, MaxFigure] = addFindFigure('CommCompos divers', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_div = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaDiversity =  taxaImpact(tSamplesNorm(:, :),HGP.Samples.topGroupName, 'diversity', samples_idx);
    tClustTaxaDiversity = tClustTaxaDiversity(tTaxaDiversity.Row, :);
    tClustTaxaDiversity.mean_diversity_dev = tClustTaxaDiversity.mean_diversity./tTaxaDiversity.mean_diversity;
    low_y = 0.35; high_y = 1.8;
    barLogWithNames(tClustTaxaDiversity.mean_diversity_dev, tTaxaDiversity.Row, clrs(k,:), [low_y,  high_y], 'YLabel','Richness deviation');
    commmun_deviations_div = [commmun_deviations_div, log10(tClustTaxaDiversity.mean_diversity_dev)];
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg,  'HShift', -0.02, 'VShift', 0.0, 'Location', 'NorthEast');
%D = (pdist(commmun_deviations,'correlation'));
% [fg, MaxFigure] = addFindFigure('Comm-similarities distributions', MaxFigure, [10, 20, 400 400], 'normal', 'on');
% histogram(D(:), 10)


%% Get most abundant phyla 
tTaxaAbundance =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'abundance', idxSamplCols);
tTaxaAbundance = tTaxaAbundance(tTaxaDiversity.Row, :);
%tTaxaAbundance = sortrows(tTaxaAbundance, "mean_diversity", "descend");
%tTaxaAbundance = tTaxaAbundance(1:HGP.Samples.topGroupNameNMostDiverse, :);

[fg, MaxFigure] = addFindFigure('CommCompos abund', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpact(tSamplesNorm(:, :),HGP.Samples.topGroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    tClustTaxaAbundance.mean_diversity_dev = tClustTaxaAbundance.mean_diversity./tTaxaAbundance.mean_diversity;
    low_y = 0.1; high_y = 10;
    barLogWithNames(tClustTaxaAbundance.mean_diversity_dev, tTaxaAbundance.Row, clrs(k,:), [low_y,  high_y], 'YLabel','Abundance deviation');
    commmun_deviations_abund = [commmun_deviations_abund, log10(tClustTaxaAbundance.mean_diversity_dev)];
    title(['Type ' NS(k)]);    
end
AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');
%D = (pdist(commmun_deviations,'correlation'));
% [fg, MaxFigure] = addFindFigure('Comm-similarities distributions', MaxFigure, [10, 20, 400 400], 'normal', 'on');
% histogram(D(:), 10)
[Lia,Locb] =ismember(tClustTaxaAbundance.Phylum, tTaxaDiversity.Phylum);

[fg, MaxFigure] = addFindFigure('CommCompos abund vs div', MaxFigure, [10, 20,400 400], 'normal', 'on');
commmun_deviations_abund = commmun_deviations_abund(Locb, :);
scatter(commmun_deviations_div(:), commmun_deviations_abund(:))
f_Lbls('Divesity deviations', 'Abundance deviations');


%% Get most abundant species 
GroupName = 'Genus';
%GroupName = 'Species';
tTaxaAbundance =  taxaImpact(tSamplesNorm,GroupName, 'abundance', idxSamplCols);
tTaxaAbundance('NA', :) = [];
tTaxaAbundance = sortrows(tTaxaAbundance, "mean_diversity", "descend");
tTaxaAbundance = tTaxaAbundance(1:20, :);
tTaxaAbundance = sortrows(tTaxaAbundance, GroupName, "ascend");


[fg, MaxFigure] = addFindFigure('CommCompos spec abund', MaxFigure, [10, 20, 1500 840], 'normal', 'on');

clrs = f_Clrs_fresh(clust_model.K);
commmun_deviations_abund = [];
for k = 1:clust_model.K
    nexttile
    idx = (tDMC.ClustID == k);
    samples_idx = idxSamplCols(idx);
    tClustTaxaAbundance =  taxaImpact(tSamplesNorm(:, :),GroupName, 'abundance', samples_idx);
    tClustTaxaAbundance = tClustTaxaAbundance(tTaxaAbundance.Row, :);
    %tClustTaxaAbundance('NA', :) = [];
    tClustTaxaAbundance.mean_diversity_dev = tClustTaxaAbundance.mean_diversity./tTaxaAbundance.mean_diversity;
    low_y = 0.05; high_y = 1/low_y;
    barLogWithNames(tClustTaxaAbundance.mean_diversity_dev, tTaxaAbundance.Row, clrs(k,:), [low_y,  high_y], 'YLabel','Abundance deviation');
    commmun_deviations_abund = [commmun_deviations_abund, log10(tClustTaxaAbundance.mean_diversity_dev)];
    title(['Clust ' NS(k)]);    
end
AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');





%% plot community type vs chlorophyll 
[fg, MaxFigure] = addFindFigure('comm_type_chl', MaxFigure, [10, 20,400 400], 'normal', 'on');
nexttile
boxplot(tStPredictors.chlorophyll_a_0, tDMC.ClustID)
ylabel('Chl-a')
xlabel('Community type')
% nexttile
% boxplot(tStPredictors.mod_aq_chl, tDMC.ClustID)
% nexttile
% ind = tStPredictors.Depth > 10;
% boxplot(tStPredictors.mod_aq_chl(ind), tDMC.ClustID(ind))
% %violinplot(tDMC.ClustID, log10(tStPredictors.chlorophyll_a_0))


% %% Get most stable species, select by smallest variance 
% % tStableSpec =  taxaImpact(tSamplesNorm,HGP.Samples.topGroupName, 'diversity', idxSamplCols);
% % tTaxaDiversity = sortrows(tTaxaDiversity, "mean_diversity", "descend");
% % tTaxaDiversity = tTaxaDiversity(1:HGP.Samples.topGroupNameNMostDiverse, :);
% 
% [fg, MaxFigure] = addFindFigure('Stable spec', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');
% 
% clrs = f_Clrs_fresh(clust_model.K);
% %dSamplesNorm is a matrix species x samples
% selectedSpecies = false(height(dSamplesNorm), 1); 
% tau = 0.8;  %minimal occupancy in the cluster 
% %remove zeros
% epsAbund = 0.5 * min(dSamplesNorm(dSamplesNorm(:)>0)); 
% dSamplesNormNoZeros = dSamplesNorm;
% dSamplesNormNoZeros(dSamplesNormNoZeros(:)==0) = epsAbund;
% 
% for k = 1:clust_model.K
%     %select a table with communities of a species in the community cluster
%     idx = (tDMC.ClustID == k);  %index of samples for cluster k
%     X = dSamplesNormNoZeros(:, idx);
%     %select specie with higher than tau occupancy 
%     selectedClusterSpecies = (sum(X(:, :)>epsAbund, 2)/width(X)) >tau;
%     selectedInd = 1:height(dSamplesNormNoZeros);
%     selectedInd = selectedInd(selectedClusterSpecies);
% 
%     %keep selected species only and normalize it again
%     X = X(selectedClusterSpecies, :);
%     X = X./ sum(X, 1);
% 
%     gm   = exp(mean(log(X), 1));   % geometric mean for each community
%     Xclr = log(X ./ gm);          % CLR (centered log-ratio) transformed abundance 
% 
%     varCLR = var(Xclr, 0, 2);                  % variance of species CLR across samples 
%     sigmaCLR = sqrt(varCLR);
%     q25 = quantile(sigmaCLR,0.05);
%     stableSpec = sigmaCLR <q25;
%     %dSamplesNorm(selectedInd(stableSpec), idx)
%     %tTraits(selectedInd(stableSpec), 'Species');
%     low_y = 1e-6; high_y = 1;
%     nexttile
%     barLogWithNames(median( dSamplesNorm(selectedInd(stableSpec), idx), 2),  tTraits.Species(selectedInd(stableSpec)), clrs(k,:), [low_y,  high_y], 'YLabel','Stable species deviation');
%     %commmun_deviations_div = [commmun_deviations_div, log10(tClustTaxaDiversity.mean_diversity_dev)];
%     title(['Type ' NS(k)]);    
% end
% AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');
% %D = (pdist(commmun_deviations,'correlation'));
% % [fg, MaxFigure] = addFindFigure('Comm-similarities distributions', MaxFigure, [10, 20, 400 400], 'normal', 'on');
% % histogram(D(:), 10)
% 
% 
% %% Get most stable species, select highest rank
% 
% [fg, MaxFigure] = addFindFigure('Stable spec', MaxFigure, [10, 20, 1.5007e+03 656], 'normal', 'on');
% 
% clrs = f_Clrs_fresh(clust_model.K);
% %dSamplesNorm is a matrix species x samples
% selectedSpecies = false(height(dSamplesNorm), 1); 
% n = 20;                     % threshold for top n species 
% tau = 0.8;  %minimal occupancy in the cluster 
% %remove zeros
% for k = 1:clust_model.K
%     %select a table with communities of a species in the community cluster
%     idx = (tDMC.ClustID == k);  %index of samples for cluster k
%     X = dSamplesNorm(:, idx);
%     [nS, nC] = size(X);
%     %select specie with higher than tau occupancy 
%     selectedClusterSpecies = (sum(X(:, :)>epsAbund, 2)/width(X)) >tau;
%     selectedInd = 1:height(dSamplesNorm);
%     selectedInd = selectedInd(selectedClusterSpecies);
% 
%     % Ranking within each community (descending), handling ties with average ranks
%     rankPos = zeros(nS, nC);            % rank 1 = most abundant
%     for j = 1:nC
%         rankPos(:, j) = tiedrank(-X(:, j));  % requires Statistics and Machine Learning Toolbox
%     end
% 
%     % Salience metrics per species
%     meanRank = mean(rankPos, 2);        % average rank across communities
%     medRank  = median(rankPos, 2);      % median rank (robust to outliers)
%     pTopK    = mean(rankPos <= K, 2);   % fraction of communities where species is in top-K
% 
%     % Optional: rank percentile (0..1), comparable across datasets; lower is better
%     rankPct  = (rankPos - 0.5) / nS;
%     meanPct  = mean(rankPct, 2);
% 
%     % Selection rules for "most salient" species
%     inTop_mean = (meanRank <= K);       % "on average in top-K"
%     inTop_freq = (pTopK >= 0.60);       % appears in top-K in ≥60% of communities
%     selected   = inTop_mean | inTop_freq;
% 
%     % Optional: Borda count within top-K (more points for higher positions)
%     bordaK = sum(max(0, K + 1 - rankPos), 2); % higher = more salient
% 
%     % Optional: assemble a results table (uses species names if available)
%     if exist('speciesNames','var')
%         rowNames = speciesNames;
%     else
%         rowNames = strcat("sp", string(1:nS));
%     end
% 
%     T = table(meanRank, medRank, pTopK, meanPct, bordaK, selected, ...
%         'RowNames', rowNames);
% 
%     % Sort by selection flag, then by Borda score (descending)
%     T = sortrows(T, {'selected','bordaK'}, {'descend','descend'});
% 
% 
% 
%     low_y = 1e-6; high_y = 1;
%     nexttile
%     barLogWithNames(median( dSamplesNorm(selectedInd(stableSpec), idx), 2),  tTraits.Species(selectedInd(stableSpec)), clrs(k,:), [low_y,  high_y], 'YLabel','Stable species deviation');
%     title(['Type ' NS(k)]);    
% end
% AddLetters2Plots(fg,  'HShift', -0.03, 'VShift', 0.0, 'Location', 'NorthEast');

