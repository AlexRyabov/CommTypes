function tTaxaDivercity = taxaImpactSamples(tSamplesNorm,topGroupName,impactType,  idxSamplCols)
%find impact values for  samples 
switch impactType
    case 'diversity'
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)>0], ...
            topGroupName, 'sum');
        tTaxaDivercity.GroupCount = [];
    case 'esn'
       esn = @(X) f_DiversityMetrics(X(:, :), [], 'SimpsonESN', [], 1);
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)], ...
           topGroupName, esn);
        tTaxaDivercity.GroupCount = [];
    case 'abundance'
        tTaxaDivercity = grpstats([tSamplesNorm(:,topGroupName), tSamplesNorm(:, idxSamplCols)], ...
            topGroupName, 'sum');
        tTaxaDivercity.GroupCount = [];
    otherwise
        error('impactType %s is undefined', impactType);
end
tTaxaDivercity.Row = tTaxaDivercity.(topGroupName);

end