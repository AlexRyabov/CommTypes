function topGroupNames = cleanSpecNames(topGroupNames)
%remove spaces and prefixes from names
topGroupNames = cellfun(@(x) strrep(x, ' ', ''), topGroupNames, 'UniformOutput', false);
% Replace 'g_' etc at the beginning of strings in the cell array
topGroupNames = cellfun(@(x) regexprep(x, '^s_', ''), topGroupNames, 'UniformOutput', false);
topGroupNames = cellfun(@(x) regexprep(x, '^g_', ''), topGroupNames, 'UniformOutput', false);
topGroupNames = cellfun(@(x) regexprep(x, '^f_', ''), topGroupNames, 'UniformOutput', false);
topGroupNames = cellfun(@(x) regexprep(x, '^o_', ''), topGroupNames, 'UniformOutput', false);
topGroupNames = cellfun(@(x) regexprep(x, '^c_', ''), topGroupNames, 'UniformOutput', false);
topGroupNames = cellfun(@(x) regexprep(x, '^p_', ''), topGroupNames, 'UniformOutput', false);


