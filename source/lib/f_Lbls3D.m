function f_Lbls3D(xLabel, yLabel, varargin)
%% Alexey Ryabov 2014
%% f_Lbls sets xlabel, ylabel and title (optional)
xlabel(xLabel);
ylabel(yLabel);
if nargin >2
    zlabel(varargin{1});
end
if nargin >3
    title(varargin{2});
end
end