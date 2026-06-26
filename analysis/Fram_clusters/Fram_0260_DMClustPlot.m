%use radial clustering and plot clusters together with diffusion map 

proj_Initialize;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
D = load(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'));
clust_model = D.clust_model;

%%


%%
[fg, fgNum] = addFindFigure('DiffMap', MaxFigure, [10, 20, 862 656], 'normal', 'on');

%az_el = [20, 80];
az_el = [0, 90];
cm = turbo;
cm = cm(end:-1:1, :);


DM_Comp2Plot = [1, 2, 3];
mark_size = 30;
t = nexttile;
scatter3(t, tDMC{:, DM_Comp2Plot(1)} , tDMC{:, DM_Comp2Plot(2)}, tDMC{:, DM_Comp2Plot(3)}, mark_size,  tStPredictors.Month_doy,  "filled", 'MarkerEdgeColor', 'k')
%view([-188.63 5.21])
colormap(t, cm)
f_Lbls3D('DM_1', 'DM_2','DM_3')
view(az_el )
colorbar(t)
title('Month')


t = nexttile;
scatter3(t, tDMC{:, DM_Comp2Plot(1)} , tDMC{:, DM_Comp2Plot(2)}, tDMC{:, DM_Comp2Plot(3)}, mark_size,  tStPredictors.Depth,  "filled", 'MarkerEdgeColor', 'k')
%view([-188.63 5.21])
colormap(t, cm)
f_Lbls3D('DM_1', 'DM_2','DM_3')
view(az_el )
colorbar(t)
title('Depth')

t = nexttile;
scatter3(t, tDMC{:, DM_Comp2Plot(1)} , tDMC{:, DM_Comp2Plot(2)}, tDMC{:, DM_Comp2Plot(3)}, mark_size,  tStPredictors.Lon,  "filled", 'MarkerEdgeColor', 'k')
%view([-188.63 5.21])
colormap(t, cm)
f_Lbls3D('DM_1', 'DM_2','DM_3')
view(az_el )
colorbar(t)
title('Longitude')


t = nexttile;
clrs = f_Clrs_fresh(clust_model.K);

for k = 1:clust_model.K
    idx = (tDMC.ClustID == k);
    s = scatter3(t, tDMC{idx,1}, tDMC{idx,2}, tDMC{idx,3}, ...
             mark_size, clrs(k,:), 'filled', 'MarkerEdgeColor', 'k');
   s.DisplayName = sprintf('C %d', k);
hold on    
end
legend('Location','bestoutside')
hold off
f_Lbls3D('DM_1', 'DM_2','DM_3')
view(az_el )
title('Month')

title('Community types')

AddLetters2Plots(fg,  'HShift', -0.05, 'VShift', -0.02)
%%