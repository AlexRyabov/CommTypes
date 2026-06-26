[fg, MaxFigure] = addFindFigure('Confusionmat', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');
nexttile
Y_train_hat = Mdl1.predict(X_train);
cm = confusionchart(Y_train, Y_train_hat, 'Normalization','row-normalized');
cm.Title = 'TRAIN. Confusion matrix (нормировка по строкам = Recall по классам)';
cm.RowSummary    = 'row-normalized';   % суммарные recall по строкам
cm.ColumnSummary = 'column-normalized';% суммарные precision по столбцам


nexttile
Y_test_hat = Mdl1.predict(X_test);
cm = confusionchart(Y_test, categorical(Y_test_hat, classes), 'Normalization','row-normalized');
cm.Title = 'TEST. Confusion matrix (нормировка по строкам = Recall по классам)';
cm.RowSummary    = 'row-normalized';   % суммарные recall по строкам
cm.ColumnSummary = 'column-normalized';% суммарные precision по столбцам

nexttile
pred_importance = MdlRF.OOBPermutedPredictorDeltaError; 
[~, idx] = sort(pred_importance, "descend");
bar(X_train.Properties.VariableNames(idx), pred_importance(idx));

disp(Mdl1.loss(X_test, Y_test));

[fg, MaxFigure] = addFindFigure('ConfusionmatYears', MaxFigure, [-1, -1,900 900], 'normal', 'on');
un_train_years = unique(Years_train);
for i = 1:length(un_train_years)
    idx_train_year = Years_train == un_train_years(i);
    nexttile
    cm = confusionchart(Y_train(idx_train_year), Y_train_hat(idx_train_year), 'Normalization','row-normalized');
    cm.Title = ['Year' NS(un_train_years(i))];
    cm.RowSummary    = 'row-normalized';   % суммарные recall по строкам
    cm.ColumnSummary = 'column-normalized';% суммарные precision по столбцам
end
