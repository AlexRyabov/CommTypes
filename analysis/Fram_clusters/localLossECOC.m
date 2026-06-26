function L = localLossECOC(Xtr, Ytr, Xte, Yte, Cost, BoxConstraint)
    tSVM = templateSVM('KernelFunction','linear','Standardize',true);
    % hopts = struct('CVPartition', cvInner, 'MaxObjectiveEvaluations', 25, ...
    %                'UseParallel', true, 'ShowPlots', false, 'Verbose', 0);

    M = fitcecoc(Xtr, Ytr, 'Learners', tSVM, 'Coding','onevsone', ...
                 'Cost', Cost, 'BoxConstraint', BoxConstraint);

    % Return classification error on the outer test fold
    L = loss(M, Xte, Yte);  % lower is better
end
