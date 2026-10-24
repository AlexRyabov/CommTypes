function hg_dm_plotDM1(DM,tStPredictors)
[fg, fgNum] = addFindFigure('DiffMap', 102, [], 'docked', 'on');
colormap(parula);
nexttile(1)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 2), DM.DMC(:, 3), 20,  tStPredictors.Month_doy,  "filled")
f_Lbls3D('DM_1', 'DM_2','DM_3')
view([0,  0, 90])

colorbar
title('Month')

nexttile(2)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 2), DM.DMC(:, 3), 20,  tStPredictors.Depth,  "filled")
colorbar
f_Lbls3D('DM_1', 'DM_2','DM_3')
view([0,  0, 90])
title('Depth')

nexttile(3)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 2), DM.DMC(:, 3), 20,  tStPredictors.Lon,  "filled")
colorbar
f_Lbls3D('DM_1', 'DM_3','DM_2')
view([0,  0, 90])
title('Lon')

nexttile(4)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 3), DM.DMC(:, 2), 20,  tStPredictors.Lat,  "filled")
colorbar
f_Lbls3D('DM_1', 'DM_3','DM_2')
view([0,  0, 90])
title('Lat')

nexttile(5)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 3), DM.DMC(:, 2), 20,  tStPredictors.temp,  "filled")
colorbar
f_Lbls3D('DM_1', 'DM_3','DM_2')
view([0,  0, 90])
title('temp')

nexttile(6)
scatter3(DM.DMC(:, 1) , DM.DMC(:, 3), DM.DMC(:, 2), 20,  tStPredictors.cond,  "filled")
colorbar
f_Lbls3D('DM_1', 'DM_3','DM_2')
view([0,  0, 90])
title('cond')
end