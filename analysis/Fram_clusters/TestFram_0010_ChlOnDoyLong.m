proj_Initialize
%%
%PlotID = 'Chl'; 
PlotID = 'T'; 

switch PlotID
    case 'Chl'
        chlFileName = [HGP.Data.Path2RawPangaeaChl, HGP.Data.fnPangaeaChl];
        field = 'chlorophyll_a_0';
        scale = 'log';
        clrLimits = [1e-3 50];
    case 'T'
        chlFileName = [HGP.Data.Path2RawPangaeaCTD,  HGP.Data.fnPangaeaCTD];
        field = 'temp';
        scale = 'linear';
        clrLimits = [-3 3];
    otherwise
        error('Undefined')
end



% Read data
T = readtable(chlFileName);

% Convert time and calculate day of year
if ~isdatetime(T.date_time)
    T.date_time = datetime(T.date_time, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
end
T.day_of_year = day(T.date_time, 'dayofyear');

% Add a small longitude jitter
rng(1);
T.long_jitter = T.long + (rand(height(T),1) - 0.5) * 0.12;

% Depth ranges
depthRanges = [0 5;
               5 15;
               15 30;
               30 60;
               60 110];

depthLabels = {'0-5 m', '5-15 m', '15-30 m', '30-60 m', '60-110 m'};

% Valid rows for plotting
% Note: chlorophyll_a_0 must be > 0 for logarithmic color scale
validRows = ~isnat(T.date_time) & ...
            ~isnan(T.long) & ...
            ~isnan(T.depth) & ...
            ~isnan(T.(field)) ;

if strcmp(scale, 'log')
    validRows = validRows & ...
            T.(field) > 0;
end 

% Common axis limits for all panels
xLimits = [min(T.long_jitter(validRows)), max(T.long_jitter(validRows))];
yLimits = [min(T.day_of_year(validRows)), max(T.day_of_year(validRows))];

% Figure
figure(1);
clf
tl = tiledlayout('TileSpacing','compact','Padding','compact');

ax = gobjects(size(depthRanges,1),1);

for i = 1:size(depthRanges,1)
    zmin = depthRanges(i,1);
    zmax = depthRanges(i,2);

    if i < size(depthRanges,1)
        rows = validRows & T.depth >= zmin & T.depth < zmax;
    else
        rows = validRows & T.depth >= zmin & T.depth <= zmax;
    end

    ax(i) = nexttile;
    
    scatter(T.long_jitter(rows), T.day_of_year(rows), ...
        35, T.(field)(rows), 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

    xlabel('Longitude');
    ylabel('Day of year');
    title(['Depth ', depthLabels{i}]);
    grid on;
    box on;

    % Logarithmic color scale
    if strcmp(scale, 'log')
    ax(i).ColorScale = 'log';
    clim(clrLimits);
    else
        ax(i).ColorScale = 'log';
        clim(clrLimits);
    end
    colormap(ax(i), jet);

    % Same x/y limits in all panels
    xlim(xLimits);
    ylim(yLimits);

    % Colorbar with log ticks
    cb = colorbar;
    cb.Ticks = [1e-3 1e-2 1e-1 1 10 50];
    cb.TickLabels = {'0.001','0.01','0.1','1','10','50'};
    cb.Label.String = 'chlorophyll a';
end

% Link axes across panels
linkaxes(ax, 'xy');

sgtitle(tl, 'Chlorophyll a by longitude and day of year');