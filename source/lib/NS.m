function s = NS(n, varargin)
%% string to num
if length(varargin) > 0
    s = num2str(n, varargin{1});
else
    s = num2str(n);
end

