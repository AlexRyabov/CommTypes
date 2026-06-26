%plot distribution of community types as a function of long and month for
%three depth layers 

Fram_0310_MLLocation_Init

%Load trained location model
tmp = load(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']), 'Mdl1');
Mdl1 = tmp.Mdl1;
clearvars tmp;



%% Create figure
[fg, MaxFigure] = addFindFigure('DominantClass', MaxFigure, [-1, -1,1.2513e+03 400], 'normal', 'on');


% depths = [0, 50; 50, 100; ];
% titles = {'Shallow (0–50 m)', 'Deep (50–100 m)'};
%depths = [0, 35; 35, 60; 60, 100];
depths = [0, 30; 30, 60; 60, 100];
titles = {'Shallow (0–30 m)', 'Middle (30-60m)', 'Deep (60–100 m)'};
Lon0 = 0;
Month0 = 7;
lonVec   = linspace(min(tStPredictors.Lon), max(tStPredictors.Lon), 50)';     %
monthVec   = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 50)';     %
[LonGrid, MonthGrid] = meshgrid(lonVec, monthVec);

fun = @(X) ml_returnPosterior(Mdl1, X);

pds = cell(height(depths));
for iD = 1:height(depths)

    depthGrid = (depths(iD, 1):5:depths(iD, 2))';
    selectedData = table( ...
        repmat(Lon0,   size(depthGrid)), ...
        repmat(Month0, size(depthGrid)), ...
        depthGrid, ...
        'VariableNames', {'Lon','Month_doy','Depth'});

    [pd,x,y] = partialDependence(fun, {'Lon', 'Month_doy'}, selectedData, ...
        "QueryPoints", [LonGrid(:) MonthGrid(:)], ...
        "UseParallel",true,  "OutputColumns", 1:8);
    pds{iD} = pd;
    K = size(pd, 1);   % number of classes (8)

    % Find max probability and corresponding class index at each grid point
    [pMax, classIdx] = max(pd, [], 1);   % size: 1 x Ny x Nx

    % Squeeze to 2D matrices: Ny x Nx
    pMax2D    = squeeze(pMax);      % max posterior at each (Month, Lon)
    classMap2D = squeeze(classIdx); % winning class index 1..K

    nexttile(iD)
    
    pcolor(x, y, classMap2D);      % x -> Lon, y -> Month_doy
    clrs = f_Clrs_fresh(K);
    colormap(clrs)
    if iD == height(depths)
        colorbar;
    end
    caxis([0.5 K+0.5]); 
    shading flat
    % color scale exactly from 1 to K
    xlabel('Longitude');
    ylabel('Month');
    title(titles(iD));
    hold on;
    contour(x, y, pMax2D, [0, 0.4, 0.6, 0.8, 1], 'k', 'LineWidth', 1, 'ShowText','on');
    hold off;
end

