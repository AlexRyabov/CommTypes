function [imp, impTable, baseLoss] = permutationImportance(X, Y, Mdl, cvp)
% X: n×p matrix (standardized or set Standardize=true in t)
% Y: labels (categorical recommended)
% cvp: cvpartition used to train CV model 
% Mdl: trained  Model

Y = categorical(Y);

% Train CV ECOC model with your partition
%CVMdl = fitcecoc(X, Y, 'Learners', t, 'CVPartition', cvp);

% Baseline CV loss (average over folds)
K = cvp.NumTestSets;
cvMdl = crossval(Mdl, 'CVPartition', cvp);
baseLoss = kfoldLoss(cvMdl); 

% Permutation importance
p = size(X,2);
imp = zeros(p,1);

for j = 1:p
    ploss = zeros(K,1);
    for k = 1:K
        Mdlk = cvMdl.Trained{k};
        idxVal = cvp.test(k);

        Xv = X(idxVal,:);
        % permute ONLY the j-th column within the validation fold
        Xv(:,j) = Xv(randperm(sum(idxVal)), j);

        ploss(k) = loss(Mdlk, Xv, Y(idxVal));
    end
    imp(j) = mean(ploss) - baseLoss;  % increase in loss ⇒ importance
end

%  summary table
impTable = table((1:p)', imp, 'VariableNames', {'PredictorIndex','PermutedLossIncrease'});
[~, order] = sort(impTable.PermutedLossIncrease, 'descend');
impTable = impTable(order,:);
end
