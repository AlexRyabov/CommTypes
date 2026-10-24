function summaryCounts = countTaxonsNumbers(tSamples)
% -----------------------------------------------------------
% Count ASVs in tSamples at each taxonomic resolution.
% Column name is **Species**.
%
% Suffix convention (PR²):
%   _sp          → genus
%   _X_sp        → family
%   _XX_sp       → order
%   _XXX_sp      → class
%   _XXXX_sp     → phylum
%
% A trailing dot may be present (e.g. "_XX_sp.").
% Rows are removed from consideration as soon as they are
% counted, so categories do not overlap.
% -----------------------------------------------------------

% Make sure the Species column is a string array
specStrings = string(tSamples.Species);

% Work mask of rows that are still uncounted
stillUnmatched = true(size(specStrings));

% Suffixes listed from longest to shortest to avoid overlap
suffixes = { '_XXXX_sp', '_XXX_sp', '_XX_sp', '_X_sp', '_sp' };
labels   = { 'Phylum',    'Class',   'Order',   'Family', 'Genus' };

% Pre-allocate count vector (+1 slot for fully resolved species)
counts = zeros(numel(labels)+1,1);
names  = [ labels, 'Species'];

% Helper function: does a string end with suffix (with or without dot)?
endsWithEither = @(str,sfx) endsWith(str,sfx) | endsWith(str, sfx + ".");

% Loop through each unresolved category
for k = 1:numel(suffixes)
    sfx  = suffixes{k};
    mask = stillUnmatched & endsWithEither(specStrings, sfx);
    counts(k)    = sum(mask);    % store count
    stillUnmatched(mask) = false;  % drop these rows from further tests
end

% Remaining rows are fully resolved to species level
counts(end) = sum(stillUnmatched);

% Build and display summary
summaryCounts = table(names', counts, ...
                      'VariableNames', {'TaxonomicRank','NumberOfASVs'});
disp(summaryCounts);
