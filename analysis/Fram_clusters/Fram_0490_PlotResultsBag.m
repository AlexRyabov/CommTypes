RowSummary    = 'total-normalized'; % 'row-normalized';   % 
ColumnSummary = 'total-normalized'; % 'column-normalized';% 


[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');
nexttile
try 
    Y_train_hat = Mdl1.predict(X_train);
catch   
    Y_train_hat = Mdl1.predict(X_train, 'BinaryLoss', 'exponential');
end
cm = confusionchart(Y_train, categorical(Y_train_hat, classes), 'Normalization','row-normalized');
cm.Title = 'TRAIN. Confusion matrix (нормировка по строкам = Recall по классам)';
cm.RowSummary    = RowSummary;   
cm.ColumnSummary = ColumnSummary;


nexttile
try
    Y_test_hat = Mdl1.predict(X_test);
catch 
    Y_test_hat = Mdl1.predict(X_test, 'BinaryLoss', 'exponential');
end
cm = confusionchart(Y_test, categorical(Y_test_hat, classes), 'Normalization','row-normalized');
cm.Title = 'TEST. Confusion matrix (нормировка по строкам = Recall по классам)';
cm.RowSummary    =  RowSummary;   
cm.ColumnSummary =  ColumnSummary;
drawnow
nexttile
 if isprop(Mdl1,'OOBPermutedPredictorDeltaError')
     pred_importance = Mdl1.OOBPermutedPredictorDeltaError;
     [~, idx] = sort(pred_importance, "descend");
     bar(cleanUndersc(X_train.Properties.VariableNames(idx)),(pred_importance(idx)));
 else
    % tImp =  permutationImportance(Mdl1);
    [pred_importance, impTable, baseLoss] = permutationImportance(X_train, Y_train, Mdl1, cv_inner);
     [~, idx] = sort(pred_importance, "descend");
 %   bar(cleanUndersc(tImp.Predictor),tImp.ImportanceMean);
      bar(cleanUndersc(X_train.Properties.VariableNames(idx)),(pred_importance(idx)));
 end


e = mean(categorical(Y_test_hat, classes) ~= Y_test); 
disp('Test error');
disp(e);

[fg, MaxFigure] = addFindFigure('ConfusionmatYears', MaxFigure, [-1, -1,900 900], 'normal', 'on');
un_train = unique(Cruise_train);
for i = 1:length(un_train)
    idx_train_year = strcmp(Cruise_train, un_train(i));
    nexttile
    cm = confusionchart(Y_train(idx_train_year), categorical( Y_train_hat(idx_train_year), classes), 'Normalization','row-normalized');
    cm.Title = ['Year' (un_train(i))];
    cm.RowSummary    =  RowSummary;   
    cm.ColumnSummary =  ColumnSummary;
end
