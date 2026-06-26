function [foldIdx, cruiseFold] = makeCruiseMonthFolds(cruise, month, K)
%MAKECRUISEMONTHFOLDS Blocked, month-diversity-aware K-fold assignment.
%
%   [foldIdx, cruiseFold] = MAKECRUISEMONTHFOLDS(cruise, month, K)
%   assigns each observation to one of K folds, under the constraints:
%       - all rows belonging to the same cruise go to the same fold
%       - each fold receives approximately the same number of cruises
%         (difference at most 1)
%       - cruises are assigned to folds to maximize the diversity of
%         months in each fold (greedy heuristic)
%
%   INPUT:
%       cruise  - vector of cruise IDs (numeric, string, cellstr,
%                 or categorical), length N (one value per observation)
%       month   - vector of month labels (numeric or categorical), length N
%       K       - number of folds (e.g. 5)
%
%   OUTPUT:
%       foldIdx    - vector of length N, fold index (1..K) for each row
%       cruiseFold - table mapping each unique cruise to its fold
%
%   NOTE:
%       This function uses a greedy heuristic:
%       1) compute the month histogram for each cruise
%       2) sort cruises by the number of distinct months (descending)
%       3) for each cruise, assign it to the fold that:
%            - still needs cruises (target cruises per fold)
%            - yields the largest increase in distinct months in that fold
%          ties are broken by choosing the fold with fewer cruises so far.

    if nargin < 3
        K = 5;
    end

    cruise = cruise(:);
    month  = month(:);

    if numel(cruise) ~= numel(month)
        error('cruise and month must have the same length.');
    end

    % Unique cruises and their indices (works for cellstr, string, numeric, categorical)
    [uCruise, ~, idxCruise] = unique(cruise, 'stable');
    nCruise = numel(uCruise);

    if nCruise < K
        error('Number of unique cruises (%d) is less than K (%d).', nCruise, K);
    end

    % Normalize months to a compact index 1..M
    [uMonth, ~, idxMonth] = unique(month, 'stable');
    nMonth = numel(uMonth);

    % Month histograms per cruise (nCruise x nMonth)
    cruiseMonthHist = zeros(nCruise, nMonth);
    for c = 1:nCruise
        rows = (idxCruise == c);
        cruiseMonthHist(c, :) = accumarray(idxMonth(rows), 1, [nMonth, 1])';
    end

    % Number of distinct months per cruise
    distinctMonthsPerCruise = sum(cruiseMonthHist > 0, 2);

    % Sort cruise indices by decreasing month diversity (harder cruises first)
    [~, sortIdx] = sort(distinctMonthsPerCruise, 'descend');
    sortedCruiseIdx = sortIdx;   % indices 1..nCruise in desired order

    % Target number of cruises per fold (balanced: difference at most 1)
    base = floor(nCruise / K);
    rem  = mod(nCruise, K);
    targetCruisesPerFold = base * ones(1, K);
    targetCruisesPerFold(1:rem) = targetCruisesPerFold(1:rem) + 1;

    % State of folds
    cruisesPerFold = zeros(1, K);
    monthHistFold  = zeros(K, nMonth);  % cumulative month hist per fold

    % Mapping cruise (unique index) -> fold
    cruiseFoldIdx = NaN(nCruise, 1);

    % Greedy assignment of each cruise (by index, no '==' on cell arrays)
    for s = 1:nCruise
        cIdx  = sortedCruiseIdx(s);  % this is an index 1..nCruise
        cHist = cruiseMonthHist(cIdx, :);

        bestFold = NaN;
        bestGain = -Inf;
        bestLoad = Inf;

        % Try all folds and pick the one with:
        %   - remaining capacity
        %   - maximum gain in distinct months
        %   - among ties, minimal number of cruises so far
        for f = 1:K
            if cruisesPerFold(f) >= targetCruisesPerFold(f)
                continue; % this fold is already at target capacity
            end

            currentHist     = monthHistFold(f, :);
            currentDistinct = sum(currentHist > 0);
            newHist         = currentHist + cHist;
            newDistinct     = sum(newHist > 0);
            gain            = newDistinct - currentDistinct;

            if gain > bestGain
                bestGain = gain;
                bestFold = f;
                bestLoad = cruisesPerFold(f);
            elseif gain == bestGain
                % tie-breaker: prefer fold with fewer cruises so far
                if cruisesPerFold(f) < bestLoad
                    bestFold = f;
                    bestLoad = cruisesPerFold(f);
                end
            end
        end

        % Fallback: if all folds are "full" by target (should rarely happen),
        % assign to the fold with minimal load.
        if isnan(bestFold)
            [~, bestFold] = min(cruisesPerFold);
        end

        % Assign this cruise index to the chosen fold
        cruiseFoldIdx(cIdx) = bestFold;
        cruisesPerFold(bestFold) = cruisesPerFold(bestFold) + 1;
        monthHistFold(bestFold, :) = monthHistFold(bestFold, :) + cHist;
    end

    % Expand from unique cruises back to all rows
    foldIdx = cruiseFoldIdx(idxCruise);

    % Also return mapping cruise -> fold as a table
    cruiseFold = table(uCruise, cruiseFoldIdx, ...
        'VariableNames', {'Cruise', 'Fold'});
end
