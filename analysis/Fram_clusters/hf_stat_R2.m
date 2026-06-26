function r2 = hf_stat_R2(observed, predicted)
%HF_STAT_R2 Coefficient of determination for observed values explained by predictions.

    idx = ~isnan(observed) & ~isnan(predicted) & isfinite(observed) & isfinite(predicted);
    observed = observed(idx);
    predicted = predicted(idx);

    if numel(observed) < 2
        r2 = NaN;
        return
    end

    ssTot = sum((observed - mean(observed, 'omitnan')).^2, 'omitnan');
    ssRes = sum((observed - predicted).^2, 'omitnan');

    if ssTot <= 0
        r2 = NaN;
    else
        r2 = 1 - ssRes / ssTot;
    end
end
