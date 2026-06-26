% Fram_0277_CoreSpecContribution_init.m
%
% Shared settings for core species chlorophyll contribution scripts.

chlorophyllTransform = "power025";

switch chlorophyllTransform
    case "power025"
        chlTransform = @(y) y .^ 0.25;
        chlBackTransform = @(y) y .^ 4;
    case "log10"
        chlTransform = @(y) log10(y);
        chlBackTransform = @(y) 10 .^ y;
    otherwise
        error('Unknown chlorophyllTransform: %s.', chlorophyllTransform);
end

modelMethod = "ridge";
%modelMethod = "stepwise";
stepwiseCriterion = "aic";
stepwisePEnter = 0.05;
stepwisePRemove = 0.10;
lambdaGrid = logspace(-2, 3, 100);
lambdaSelectionTolerance = 0.05;
nFolds = 5;
blockedCV = true;
holdoutFraction = 0.20;
useClusterBalancedRidgeWeights = true;
contributionCILevel = 0.5;
contributionQuantileProbs = [0.25 0.50 0.75];
contributionIntervalLabel = "Q25Q75";
rng(1);

% Options: "none", "clr", "power025", or "power05".
modelAbundanceTransform = "power025";

roleClrs = containers.Map;
roleClrs('Autotroph')   = [0.08 0.46 0.20];
roleClrs('Mixotroph')   = [0.36 0.25 0.70];
roleClrs('Zooplankton') = [0.82 0.38 0.08];
roleClrs('Parasite')    = [0.78 0.05 0.45];
roleClrs('Recycler')    = [0.20 0.36 0.46];
roleClrs('Unknown')     = [0.15 0.15 0.15];

roleShortNames = containers.Map;
roleShortNames('Autotroph')   = 'Auto';
roleShortNames('Mixotroph')   = 'Mixo';
roleShortNames('Zooplankton') = 'Zoo';
roleShortNames('Parasite')    = 'Para';
roleShortNames('Recycler')    = 'Recy';
roleShortNames('Unknown')     = 'Unk';
