%Test effects of the threshold on the min mean abundance on diffusion map 
% and on the predictability of clusters for a chosen ML model 
clearvars
proj_Initialize;

%% Load data
minOccurence =  logspace(-2, -0.5, 20);
minMedianInPres =  [1e-5, 1e-4, 1e-3, 1e-2];
vals = minMedianInPres;
test_error = [];
for iMO = 1:length(vals)
    HGP = projParams(JobID, DataParams);
    HGP.Samples.minOccurence = 0.15; %0.024;
    HGP.Samples.minMedianInPresence = vals(iMO);
    [Data, DI] = dt_loadData(HGP, []);
    tSamplesNorm = Data.tSamplesNorm;
    tStPredictors = Data.tStPredictors;

    %%

    %% DM

    %%find the best DM parameters for which we get the largest correlation with Depth and doy
    %DM = hg_cd_aux_FindBestDMParamters(k_min,Laplacian,sigma,alpha,DMvariables,SpeciesAbund,tStPredictors);
    %metrics spearman, mcc, jaccard
    %Laplacian='rownorm'
    %k_min=55
    k_min=60;
    sigma=1.3;
    alpha = 1;
    Laplacian="rownorm";
    %Metric="jaccard";
    Metric="spearman";
    DMvariables=15;

    mask = [];

    % Find the histogram of average species abundance
    SpeciesAbundTransp = array2table(Data.dSamplesNorm');
    spNames = cleanSpecNames(tSamplesNorm.TopName);
    SpeciesAbundTransp.Properties.VariableNames = spNames;
    options = diffmap_options('k_min',k_min, 'Laplacian',Laplacian,  'Metric',Metric, 'sigma',sigma, 'alpha',alpha, 'DMvariables',DMvariables, 'mask', mask);
    DM = diffmap(SpeciesAbundTransp{:, :}, options);

    tDMC = array2table(DM.DMC);
    variableNames = arrayfun(@(x) sprintf('DM%d', x), 1:width(tDMC), 'UniformOutput', false);
    tDMC.Properties.VariableNames = variableNames;
    mkdir(HGP.Data.Path2ResultsData);
    writenewtable(tDMC, fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC.xlsx'));

    Fram_0250_Cluster;
    % Fram_0400_Init;
    % t = templateTree('MaxNumSplits',2, 'MinLeafSize', 10);  %adding  'Surrogate','all' makes it worse 
    % Mdl1 = fitcensemble(X_train, Y_train,'Learners',t, 'Method','AdaBoostM2',  'LearnRate', 0.063, ...
    % 'NumLearningCycles', 39);
    % cvMdl1= crossval(Mdl1, 'CVPartition', cv_inner);
    % test_error(iMO) =   kfoldLoss(cvMdl1);
    % 
    % [fg, MaxFigure] = addFindFigure('ErrorOnMicOcc', MaxFigure, [-1, -1,400 400], 'normal', 'on');
    % semilogx(vals(1:iMO), test_error);
    % drawnow

end

