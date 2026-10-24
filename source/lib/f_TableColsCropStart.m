function T = f_TableColsCropStart(T, startStr)
%crop the substring startStr at the beginning of column name 
%if startStr is 'mean_', then rename table column 'mean_ValueX' to'ValueX'
lMax = length(startStr);
for iC = 1:length(T.Properties.VariableNames)
    s = T.Properties.VariableNames{iC};
    if length(s) >lMax
        sStartWith = s(1:lMax);
        if strcmp(sStartWith, startStr)
             T.Properties.VariableNames{iC} = s(lMax+1:end);
        end
    end
end