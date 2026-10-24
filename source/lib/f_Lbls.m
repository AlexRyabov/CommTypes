function f_Lbls(xLabel, yLabel, varargin)
%% Alexey Ryabov 2014
%% f_Lbls sets xlabel, ylabel and title (optional)
xlabel(xLabel);
ylabel(yLabel);
if ~isempty(varargin)
    title(varargin{1})
end