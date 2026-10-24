function f_FramePlot(varargin)
%f_FramePlot('out')
%f_FramePlot('in')
box on
switch nargin 
    case 0
        set(gca, 'YMinorTick','on','XMinorTick','on','TickDir','in', 'Layer','top', 'TickLength',[0.02, 0.02]);
    case 1
        set(gca, 'YMinorTick','on','XMinorTick','on','TickDir',char(varargin{1}), 'Layer','top', 'TickLength',[0.02, 0.02]);
    case 2
        set(gca, 'YMinorTick','on','XMinorTick','on','TickDir',char(varargin{1}), 'Layer','top', 'TickLength',[varargin{2}, varargin{2}]);
end
