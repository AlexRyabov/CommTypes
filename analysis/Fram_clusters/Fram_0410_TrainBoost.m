

Fram_040_Init

%%
if 1  %tune automatically 
   % t = templateTree('MaxNumSplits',1, 'MinLeafSize', 1);
    AutoModel =  fitcensemble( ...
        X_train, Y_train, ...
          'OptimizeHyperparameters', {'NumLearningCycles','LearnRate', 'MaxNumSplits', 'MinLeafSize'}, ...
        'HyperparameterOptimizationOptions', struct( ...
        'CVPartition', cv_inner, ...
        'MaxObjectiveEvaluations', 100, ...           %
        'AcquisitionFunctionName','expected-improvement-plus', ...
        'UseParallel', true, ...
        'Verbose', 1));
     % 'Learners',t,...
        %'Method','LSBoost', ...

    disp("Auto model training " + AutoModel.loss(X_train, Y_train));
    disp("Auto model test " + AutoModel.loss(X_test, Y_test));

    % random forest 
    %    t = templateTree('MaxNumSplits',1, 'MinLeafSize', 1);
    % AutoModel =  fitrensemble( ...
    %     X_train, Y_train, ...
    %     'Learners',t,...
    %     'Method','Bag', ...
    %     'OptimizeHyperparameters', {'NumLearningCycles','MaxNumSplits', 'NumVariablesToSample'}, ...
    %     'HyperparameterOptimizationOptions', struct( ...
    %     'CVPartition', cv_inner, ...
    %     'MaxObjectiveEvaluations', 100, ...           %
    %     'AcquisitionFunctionName','expected-improvement-plus', ...
    %     'UseParallel', true, ...
    %     'Verbose', 1));
    % 
    % disp("Auto model training " + AutoModel.loss(X_train, Y_train));
    % disp("Auto model test " + AutoModel.loss(X_test, Y_test));
end


%% hyperparameter tuning
optLeaningRate = 0.1;
optNumLearningCycles = 60;
%t = templateTree('MaxNumSplits',1, 'Surrogate','all');
%t = templateTree('MaxNumSplits',1);
t = templateTree('MaxNumSplits',2, 'MinLeafSize', 10);  %adding  'Surrogate','all' makes it worse 
if 1 %%set to 1 to perform the parameter tuning
    learnRates = logspace(-3, 0, 6);
    if ~ismember(optLeaningRate, learnRates)
        learnRates = sort([learnRates, optLeaningRate]);
    end
    idx_LR = find(learnRates == optLeaningRate, 1);
    crls = f_Clrs_fresh(length(learnRates));
    [fg, MaxFigure] = addFindFigure('CVerror', MaxFigure, [-1, -1,400 400], 'normal', 'on');
    clf
    for iLR = 1:length(learnRates)
        %Mdl = fitcensemble(X_train, Y_train,'Learners',t, 'Method', 'LPBoost', 'CVPartition', cv_inner,  ...
        %    'NumLearningCycles', 1000);
         Mdl = fitcensemble(X_train, Y_train,'Learners',t, 'Method', 'AdaBoostM2', 'CVPartition', cv_inner, 'LearnRate', learnRates(iLR), ...
             'NumLearningCycles', 500);
        kflc = kfoldLoss(Mdl,'Mode','cumulative');
        figure(fg);
        learCycles = 1:length(kflc);
        p = semilogx(learCycles, kflc, 'LineWidth', 2, 'Color', crls(iLR, :));
        p.DisplayName = ['\alpha=' num2str(learnRates(iLR), 2)];
        hold on
        if idx_LR == iLR
            s = scatter(optNumLearningCycles, kflc(optNumLearningCycles), 'filled');
            s.DisplayName =['Suggested params'];
        end
        legend("Location","best")
        drawnow
        ylabel('cross-validated loss');
        xlabel('Num of learning cycles (trees) ');
    end
    
    %ylim([0.1, 1])
 

end
%%
% Mdl1 = fitcensemble(X_train, Y_train,'Learners',t, 'Method','AdaBoostM2',  'LearnRate', 0.063, ...
%     'NumLearningCycles', 39);
Mdl1 = fitcensemble(X, Y,'Learners',t, 'Method','AdaBoostM2',  'LearnRate', 1, ...
    'NumLearningCycles', 1000);


% Mdl1 = AutoModel;
Fram049_PlotResults