function tOut = cleandata_unifyNames(T, name2Unify, commonName, newName)
arguments
    T 
    name2Unify 
    commonName 
    newName = commonName; 
end
% Replace all name2Unify entries that start with commonName by newName

tOut = T;

% Logical index: true if GroupName starts with commonName
idx = startsWith(T.(name2Unify), commonName);

% Replace matched entries
tOut.(name2Unify)(idx) = newName;
end