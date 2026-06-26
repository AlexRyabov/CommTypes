%% Chl distribution as a function of Lon, Month, and Depth and Temp
%Check if adding temperature allows us to remove the dependence on the
%community type
%Conclusion 
%if we add temperature, then the dependence on the type becomes even cleaner 
%type 2 is underestimated and types 1, 5, 7, 8 are over estimated 

% Alex Ryabov, 2025

proj_Initialize

predictors = {'Lon', 'Month_doy', 'Depth', 'temp'};

%Load diffusion map clusters 
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC_clustID.xlsx'));
tStPredictors.ClustID = categorical(tDMC.ClustID);



%Remove nans
tStPredictors = tStPredictors(~isnan(tStPredictors.chlorophyll_a_0), :);


Y = tStPredictors.chlorophyll_a_0;
X = tStPredictors(:, predictors);
%% transform the response
alpha =0.25;
%alpha =1;
var_tranf = @(y) y.^alpha ;
inv_var_tranf = @(tranf_y) tranf_y.^(1/alpha);
%[Y, lambda] = boxcox(tStPredictors.chlorophyll_a_0);
% var_tranf = @(y) log10(y);
%inv_var_tranf = @(tranf_y) 10.^tranf_y;
Y = var_tranf(tStPredictors.chlorophyll_a_0);
figure(1)
clf
histogram(Y)
%%

%train svm with automatic feature selection
% [X_train, Y_train, X_test, Y_test, cv_inner] = ...
%     ml_blockedCV(X,Y, tStPredictors.Cruise, tStPredictors.MonthNum, HGP.ML.KFolds);

cv_inner = cvpartition(tStPredictors.MonthNum, "Kfold", HGP.ML.KFolds, 'Stratify',true);
nFolds = cv_inner.NumTestSets;

%Test and training
X_train = X(:, :);
Y_train = Y(:);
X_test = [];
Y_test = []; 

%%
LM_models = {'RF', 'SVM_auto', 'SVM'}'; %, 'Boost1', 'Boost2',  'Lin_stepwise', 'GAM', 'FitNet', 'RFNative'}';
Select_SaveBest_Model = false;

if Select_SaveBest_Model
    %%tune the a regressor
    Models = table(LM_models);
    Models.Row = LM_models;
    %Models.LM_models = LM_models;
    for iM = 1:length(LM_models)
        iM
        [Mdl, TuneRes] = ml_tune_model(LM_models{iM},X_train,Y_train,cv_inner);
        Models.Mdl{iM} = Mdl;
        Models.TuneRes{iM} = TuneRes;
        Models.cvloss(iM) = TuneRes.bestCVLoss;
        %plot partial dependence 
        [fg, MaxFigure] = addFindFigure(['PD' LM_models{iM}], MaxFigure, [10, 20,1200 400], 'normal', 'on');
        for iP = 1:length(predictors)
            nexttile
            plotPartialDependence( Models.Mdl{iM}, predictors{iP});
        end
    end
    save(fullfile(HGP.Data.Path2ResultsData, ['trainedChlRegrModelsLMDTemp.mat']), 'Models');
else
    tmp = load(fullfile(HGP.Data.Path2ResultsData, ['trainedChlRegrModelsLMDTemp.mat']));
    Models = tmp.Models;
    clearvars tmp;
    BestModel = Models.Mdl{'SVM_auto'};
end
   

   %plot partial dependence 
        % [fg, MaxFigure] = addFindFigure(['PD' Models.LM_models{iM}], MaxFigure, [10, 20,1200 400], 'normal', 'on');
        % for iP = 1:length(predictors)
        %     nexttile
        %     plotPartialDependence( Models.Mdl{iM}, predictors{iP});
        % end




%% plot predicted vs observed data 
[fg, MaxFigure] = addFindFigure('Chl_data_res', MaxFigure, [-1, -1, 1261 377], 'normal', 'on');


chlObs = (tStPredictors.chlorophyll_a_0);

chlHat = inv_var_tranf(BestModel.predict(tStPredictors));


% Plot predicted vs observed data
t = nexttile
%colormap(t, f_Clrs_fresh(HGP.dmCluster.NumComponents));
scatter(chlObs, chlHat, 30,  -tStPredictors.Depth, 'filled');
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
hold on;
plot([min(chlObs), max(chlObs)], [min(chlObs), max(chlObs)], 'r--'); % 45-degree line
hold off;
xlabel('Observed chlorophyll-a ');
ylabel('Predicted chlorophyll-a ');
f_Set_scales({'logx', 'logy'})
%axis equal;
%grid on;
axis([0, 5, 0, 5])
cb = colorbar
cb.Label.String = 'Depth (m)'; 

t2 = nexttile(2)
residuals = var_tranf(tStPredictors.chlorophyll_a_0) -(chlHat)  ;
scatter(tStPredictors.chlorophyll_a_0, residuals, 30,  -tStPredictors.Depth, 'filled');
xlabel('Observed chlorophyll-a (transformed)');
ylabel('Residual (Obs -Pred)');
%colormap(t2, "jet")
%caxis([0.5, HGP.dmCluster.NumComponents+0.5])
cb = colorbar
cb.Label.String = 'Depth (m)'; 
%%
% 
% 

% Calculate the residuals


%residuals = (tStPredictors.chlorophyll_a_0) - inv_var_tranf(chlHat);
residuals = var_tranf(tStPredictors.chlorophyll_a_0) - (chlHat);
% Pairwise comparisons (Tukey-Kramer)

nexttile(3)
[p, tblA, stats] = anova1(residuals, tStPredictors.ClustID, 'off');  % 'on' shows built-in plot
c = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
% c = [g1 g2 lowerCI meanDiff upperCI pValue]
disp(p);
cMulty = multcompare(stats, 'CType','tukey-kramer', "Alpha" ,0.05 );
sig = cMulty(cMulty(:,6) < 0.05, :)


p_levene = vartestn(residuals, tStPredictors.ClustID, 'TestType','LeveneAbsolute', 'Display','off');
p_levene



p_kw = kruskalwallis(residuals, categorical(tStPredictors.ClustID), 'off');
p_kw

Type = categorical(tStPredictors.ClustID);
tbl = table(residuals, Type, 'VariableNames', {'Resid','Type'});
lm = fitlm(tbl, 'Resid ~ Type', 'RobustOpts','on');
anova(lm,'summary')



