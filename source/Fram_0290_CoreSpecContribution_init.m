% Fram_0290_CoreSpecContribution_init.m
%
% Shared settings for core species chlorophyll contribution scripts.

chlorophyllTransform = "power025";

switch chlorophyllTransform
    case "none"
        chlTransform = @(y) y ;
        chlBackTransform = @(y) y ;
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
stepwiseCriterion = "aic";  %"aic" ,"BIC" "SSE"
 stepwisePEnter = 0.05;
 stepwisePRemove = 0.10;
lambdaGrid = logspace(-2, 1, 100);
lambdaSelectionTolerance = 0.0;
nFolds = 6;
blockedCV = false;
holdoutFraction = 0.30;
%There are few samples with extremely low chl values, but their exclusion does not help 
excludeLowChlorophyllTypes = false;
lowChlorophyllExcludedTypes = 2:4;
lowChlorophyllExclusionThreshold = 1e-1;
useClusterBalancedRidgeWeights = true;
contributionCILevel = 0.5;
contributionQuantileProbs = [0.25 0.50 0.75];
contributionIntervalLabel = "Q25Q75";
rng(42);

% Options: "none", "log10", "power025", or "power05".
modelAbundanceTransform = "power025";
modelAbundanceLogPseudocount = 1e-6;

[roleClrs, roleShortNames] = Fram_roleColorsAndNames();
