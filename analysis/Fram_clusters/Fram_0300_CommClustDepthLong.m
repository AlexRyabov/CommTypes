%

proj_Initialize;
% 
% Fram_010_DM
% Fram_025_Clust;

%% Load data

tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
D = load(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'));
clust_model = D.clust_model;

%%


%%
%rng(42)% to have the same random jitter
[fg, MaxFigure] = addFindFigure('Comm_distributions', MaxFigure, [10, 20, 862 656], 'normal', 'on');
plotXYByPeriods([6:10], {'June','July','August','September'}, tStPredictors.Month_doy, tStPredictors{:, {'Lon', 'Depth'}}, tDMC.ClustID)
AddLetters2Plots( 'HShift', -0.05, 'VShift', -0.02)


[fg, MaxFigure] = addFindFigure('Comm_distributions_june', MaxFigure, [10, 20,1.4687e+03 303.3333], 'normal', 'on');
plotXYByPeriods([6:1/4:7], {'June, week 1','June, week 2','June, week 3', 'June, week 4'}, tStPredictors.Month_doy, tStPredictors{:, {'Lon', 'Depth'}}, tDMC.ClustID)


[fg, MaxFigure] = addFindFigure('Comm_distributions_july', MaxFigure, [10, 20,1.4687e+03 303.3333], 'normal', 'on');
plotXYByPeriods([7:1/4:8], {'July, week 1','July, week 2','July, week 3', 'July, week 4'}, tStPredictors.Month_doy, tStPredictors{:, {'Lon', 'Depth'}}, tDMC.ClustID)

