%% Ice distribution in different months in Hausgarten
% Alex Ryabov, 2025

%Clean variables, close windows, initialize all project constants and parameters
 proj_Initialize;
% 
% 
% %% Load data
% [Data, DI] = dt_loadData(HGP, []);
% tSamplesNorm = Data.tSamplesNorm;
% tStPredictors = Data.tStPredictors;

%%
%% plot ice data
tICE = readtable(HGP.Data.fnIceMonthlyLotLat);
%it should be a table as 
% tICE =
% 
%   345870×8 table
% 
%        date        long       lat      sea_ice              ID              Month    Year    MonthName
%     __________    _______    ______    _______    ______________________    _____    ____    _________
% 
%     2009-05-15    -1.9415    78.001    27.123     {'awi.sea_ice_portal'}      5      2009     May     
%     2009-05-15     1.9415    78.001         0     {'awi.sea_ice_portal'}      5      2009     May     
%     2009-05-15     10.739    78.001         0     {'awi.sea_ice_portal'}      5      2009     May     

%%1) Make sure 'date' is datetime and add Month fields
if ~isdatetime(tICE.date)
    % adjust InputFormat if your date strings differ
    tICE.date = datetime(tICE.date,'InputFormat','yyyy-MM-dd');
end
tICE.Month     = month(tICE.date);                 % numeric month 1..12
tICE.Year     = year(tICE.date);                 % numeric month 1..12
% year_from = 2009;
% year_to = 2021;
% tICE = tICE(tICE.Year>=year_from & tICE.Year<2021, :);
tICE.MonthName = categorical(month(tICE.date,'name'));  % e.g., "June"

tStList = readtable(HGP.Data.fnStatCoors);

%%2) Average sea ice by Month, long, lat using grpstats
G = grpstats(tICE, {'Month','long','lat'}, 'mean', 'DataVars','sea_ice');
% grpstats names the averaged variable 'mean_sea_ice'
G.Properties.VariableNames{'mean_sea_ice'} = 'ice_mean';

%%3) Plot June–September in a tiled layout
monthsToPlot = 6:9;  % Jun, Jul, Aug, Sep
monthNames   = month(datetime(2000,monthsToPlot,1),'name'); % for titles

% Region limits to zoom to your data
latlim = [78,  81];
lonlim = [min(G.long)-1, max(G.long)+1];

% Custom colormap: deep-sea blue (0%) -> white (100% ice)
nC = 256;
sea = [0 0.2 0.6]; ice = [1 1 1];
t = linspace(0,1,nC)';                 % interpolation parameter
cmap = sea.*(1-t) + ice.*t;            % gradient from blue to white

MaxFigure = f_MakeFigure(MaxFigure, [-1, -1, 867.3333 668]);
orange = [1, 0.498, 0];
%boundaries of ice edge 
tIceEdgeLoc = table;

tL = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for k = 1:numel(monthsToPlot)
    m  = monthsToPlot(k);
    Tk = G(G.Month==m, :);
    ax = geoaxes(tL);              % create a geographic axes in the layout
    ax.Layout.Tile = k;           % send it to tile k
    %geodensityplot(Tk.lat, Tk.long,Tk.ice_mean)
    geoscatter(ax, Tk.lat, Tk.long, 18, Tk.ice_mean, 'filled', ...
        'MarkerEdgeColor','none');
    geobasemap('bluegreen');        % clean land/sea background
    geolimits(ax, latlim, lonlim);


    % --- Build a grid from scattered points (lon/lat -> grid) ---
    lonv = linspace(min(Tk.long), max(Tk.long), 200);
    latv = linspace(min(Tk.lat),  max(Tk.lat),  201);
    [Lon,Lat] = meshgrid(lonv, latv);
    IceGrid = griddata(Tk.long, Tk.lat, Tk.ice_mean, Lon, Lat, 'natural'); % or 'linear'

    % --- Compute the 20% contour without plotting ---
    C = contourc(lonv, latv, IceGrid, [20 20]);  % returns contour matrix

    % --- Parse contour segments and plot them on geoaxes ---
    hold on
    ii = 1;
    while ii < size(C,2)
        n   = C(2,ii);              % number of points in this segment
        seg = C(:, ii+1:ii+n);      % [lon; lat]
        lon = seg(1,:); lat = seg(2,:);
        geoplot(ax, lat, lon, 'r-', 'LineWidth', 1.8); % red 15% line
        ii = ii + n + 1;
        if n > 100  %is is a large segment 
            IceEdgeLoc.Month = m;
            ind =  78.7 < lat & lat < 79.8;
            IceEdgeLoc.ICE20LongMin = min(lon(ind));
            IceEdgeLoc.ICE20LongMax = max(lon(ind));
            tIceEdgeLoc = [tIceEdgeLoc; struct2table(IceEdgeLoc)];
        end 
    end

    geoscatter(tStList.mean_lat_dec, tStList.mean_lon_dec, 20, orange, 'filled')

    hold off

    title(ax, char(monthNames(k)));

    colormap(ax, cmap);
    caxis(ax, [0 100]);                % 0..100 % ice concentration
    cb = colorbar(ax);
    cb.Label.String = 'Sea ice concentration (%)';
end
sgtitle('Monthly mean sea ice (June–September)');

writenewtable(tIceEdgeLoc,  fullfile(HGP.Data.Path2ResultsData, 'Ice_Edge.xlsx'))
