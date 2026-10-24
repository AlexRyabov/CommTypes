% Plot distribution of community types as a function of longitude and month
% for three depth layers.
%
% calcPlot = true  -> calculate partial dependencies using the trained model
%                    and save results
% calcPlot = false -> load previously saved partial dependences and only plot

Fram_0310_SVM_Type_Init_LonMonth;
calcPlot = true;
% This script uses the LonMonthDepth predictor set.
% For plotting PD in Temperature-DOY axes use Fram_0354_ML_Plot_TempDepDoy

%% Parameters

depths = [0, 30; 
          30, 60; 
          60, 100];

titles = {'Shallow (0–30 m)', ...
          'Middle (30–60 m)', ...
          'Deep (60–100 m)'};

Lon0   = 0;
Month0 = 7;

nGrid = 20;

% File names for loading the trained model and saving/loading PD results
modelFileName = fullfile( ...
    HGP.Data.Path2ResultsData, ...
    ['trainedLocationModel_' LM_Model '.mat']);

pdFileName = fullfile( ...
    HGP.Data.Path2ResultsData, ...
    ['PD_ClustIDonLonDoyDepth_' LM_Model '.mat']);


   tmp = load(modelFileName, 'CVStats');
    CVStats = tmp.CVStats;
disp( CVStats.BestParams)

%% Calculate or load partial dependences

if calcPlot

    % Load trained location model
    tmp = load(modelFileName, 'Mdl1');
    Mdl1 = tmp.Mdl1;
    clearvars tmp;

    % Define grid
    lonVec = linspace( ...
        min(tStPredictors.Lon), ...
        max(tStPredictors.Lon), ...
        nGrid)';

    monthVec = linspace( ...
        min(tStPredictors.Month_doy), ...
        max(tStPredictors.Month_doy), ...
        nGrid)';

    [lonGrid, monthGrid] = meshgrid(lonVec, monthVec);

    % Function returning class posterior probabilities
    fun = @(X) ml_returnPosterior(Mdl1, X);

    % Storage
    pds = cell(size(depths, 1), 1);
    xs  = cell(size(depths, 1), 1);
    ys  = cell(size(depths, 1), 1);

    for iD = 1:size(depths, 1)

        depthGrid = (depths(iD, 1):5:depths(iD, 2))';

        selectedData = table( ...
            repmat(Lon0,   size(depthGrid)), ...
            repmat(Month0, size(depthGrid)), ...
            depthGrid, ...
            'VariableNames', {'Lon', 'Month_doy', 'Depth'});

        [pd, x, y] = partialDependence( ...
            fun, ...
            {'Lon', 'Month_doy'}, ...
            selectedData, ...
            "QueryPoints", [lonGrid(:), monthGrid(:)], ...
            "UseParallel", true, ...
            "OutputColumns", 1:8);

        pds{iD} = pd;
        xs{iD}  = x;
        ys{iD}  = y;
    end

    % Save everything needed for plotting
    save(pdFileName, ...
        'pds', 'xs', 'ys', ...
        'depths', 'titles', ...
        'Lon0', 'Month0', ...
        'lonVec', 'monthVec', ...
        'nGrid');

else

    % Load previously calculated partial dependencies
    tmp = load(pdFileName, ...
        'pds', 'xs', 'ys', ...
        'depths', 'titles');

    pds    = tmp.pds;
    xs     = tmp.xs;
    ys     = tmp.ys;
    depths = tmp.depths;
    titles = tmp.titles;

    clearvars tmp;

end

%% Create figure

[fg, MaxFigure] = addFindFigure( ...
    'DominantClass', ...
    MaxFigure, ...
    [-1, -1, 1.2513e+03, 400], ...
    'normal', ...
    'on');

tiledlayout(1, size(depths, 1), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for iD = 1:size(depths, 1)

    pd = pds{iD};
    x  = xs{iD};
    y  = ys{iD};

    K = size(pd, 1);   % number of classes

    % Find maximum posterior probability and winning class
    [pMax, classIdx] = max(pd, [], 1);

    % Convert to 2D matrices
    pMax2D     = squeeze(pMax);
    classMap2D = squeeze(classIdx);

    nexttile(iD)

    pcolor(x, y, classMap2D);      % x -> Lon, y -> Month_doy

    clrs = f_Clrs_fresh(K);
    colormap(clrs)

    if iD == size(depths, 1)
        colorbar;
    end

    caxis([0.5, K + 0.5]);
    shading flat

    xlabel('Longitude');
    ylabel('Month');
    title(titles{iD});

    hold on;

    % Plot grey transparent mask outside the observed data cloud
    indDepth = tStPredictors.Depth >= depths(iD, 1) & ...
               tStPredictors.Depth <  depths(iD, 2);

    selectedPreds = tStPredictors(indDepth, :);

    xData = selectedPreds.Lon;
    yData = selectedPreds.Month_doy;
    clustData = tDMC.ClustID(indDepth);

    % Remove NaNs
    indGood = ~isnan(xData) & ~isnan(yData);
    xData = xData(indGood);
    yData = yData(indGood);
    clustData = clustData(indGood);

    xLimits = [min(x(:)), max(x(:))];
    yLimits = [min(y(:)), max(y(:))];

    xlim(xLimits);
    ylim(yLimits);

    if numel(xData) >= 3

        % Outer rectangle covering the whole plotting area
        outerPoly = polyshape( ...
            [xLimits(1), xLimits(2), xLimits(2), xLimits(1)], ...
            [yLimits(1), yLimits(1), yLimits(2), yLimits(2)]);

        % Hull around observed data
        shrinkFactor = 0.2;   % 1 = convex hull, smaller = more concave boundary
        k = boundary(xData, yData, shrinkFactor);

        dataPoly = polyshape(xData(k), yData(k));

        % Area outside observed data
        maskPoly = subtract(outerPoly, dataPoly);

        plot(maskPoly, ...
            'FaceColor', [0.7, 0.7, 0.7], ...
            'FaceAlpha', 0.5, ...
            'EdgeColor', 'none');
    end

    % Probability contours
    contour(x, y, pMax2D, ...
        [0, 0.4, 0.6, 0.8, 1], ...
        'k', ...
        'LineWidth', 1, ...
        'ShowText', 'on');

    % Plot original observations last so they stay visible
    scatter(xData, yData, ...
        20, ...
        clustData, ...
        'filled', ...
        'MarkerEdgeColor', 'k');

    hold off;

end
AddLetters2Plots(fg, {'a', 'b', 'c'},  HShift=-0.02, VShift=-0.1)
