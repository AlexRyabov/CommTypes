function [tTaxaDivercity, tTaxaDivercitySamples] = taxaImpact(tSamplesNorm,topGroupName,impactType,  idxSamplCols)
%find median across samples
switch impactType
    case 'diversity'
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)>0], ...
            topGroupName, 'sum');
        tTaxaDivercity.GroupCount = [];
        tTaxaDivercitySamples = tTaxaDivercity;
        tTaxaDivercitySamples.Row =  tTaxaDivercitySamples.(topGroupName);
        tTaxaDivercitySamples.(topGroupName) = [];
        tTaxaDivercity  = table(tTaxaDivercity.(topGroupName), median(tTaxaDivercity{:, 2:end}, 2), quantile(tTaxaDivercity{:, 2:end}, [0.25], 2),  quantile(tTaxaDivercity{:, 2:end}, [0.75], 2), std(tTaxaDivercity{:, 2:end},0, 2));
        tTaxaDivercity.Properties.VariableNames = {topGroupName, 'mean_diversity', 'quantile25', 'quantile75', 'std'};
    case 'esn'
        esn = @(X) f_DiversityMetrics(X(:, :), [], 'SimpsonESN', [], 1);
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)], ...
            topGroupName, esn);
        tTaxaDivercity.GroupCount = [];
        tTaxaDivercitySamples = tTaxaDivercity;
        tTaxaDivercitySamples.Row =  tTaxaDivercitySamples.(topGroupName);
        tTaxaDivercitySamples.(topGroupName) = [];
        tTaxaDivercity  = table(tTaxaDivercity.(topGroupName), median(tTaxaDivercity{:, 2:end}, 2), quantile(tTaxaDivercity{:, 2:end}, 0.25, 2),  quantile(tTaxaDivercity{:, 2:end}, 0.75, 2), std(tTaxaDivercity{:, 2:end},0, 2));
        tTaxaDivercity.Properties.VariableNames = {topGroupName, 'mean_diversity', 'quantile25', 'quantile75', 'std'};
    case 'abundance'
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)], ...
            topGroupName, 'sum');
        tTaxaDivercity.GroupCount = [];
        tTaxaDivercitySamples = tTaxaDivercity;
        tTaxaDivercitySamples.Row =  tTaxaDivercitySamples.(topGroupName);
        tTaxaDivercitySamples.(topGroupName) = [];
        tTaxaDivercity  = table(tTaxaDivercity.(topGroupName), median(tTaxaDivercity{:, 2:end}, 2), quantile(tTaxaDivercity{:, 2:end}, 0.25, 2),  quantile(tTaxaDivercity{:, 2:end}, 0.75, 2), std(tTaxaDivercity{:, 2:end},0, 2));
        tTaxaDivercity.Properties.VariableNames = {topGroupName, 'mean_diversity', 'quantile25', 'quantile75', 'std'};
    otherwise
        error('impactType %s is undefined', impactType);
end
tTaxaDivercity.Row = tTaxaDivercity.(topGroupName);

end