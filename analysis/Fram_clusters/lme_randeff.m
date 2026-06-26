function [selStat, selStd] = lme_randeff(mle, coefName, varName)
%Get statistics and std for given coefficient name and group of random
%effects
[~, ~, Stat] = randomEffects(mle, 'Alpha', 0.05); 
isInt   = strcmp(Stat.Name, coefName) & strcmp(Stat.Group, varName);
selStat = Stat(isInt, :);
selStd = std(selStat.Estimate);

end