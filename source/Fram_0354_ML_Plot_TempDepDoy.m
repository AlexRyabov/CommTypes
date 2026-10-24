% Plot distribution of community types as a function of temperature and month
% for three depth layers.
%
Fram_0350_ML_Init_TempDepDoy

calcPlot = true;
% calcPlot = true  -> calculate partial dependences using the trained model
%                    and save results
% calcPlot = false -> load previously saved partial dependences and only plot



%% Parameters

depths = [0, 30; 
          30, 60; 
          60, 100];

titles = {'Shallow (0–30 m)', ...
          'Middle (30–60 m)', ...
          'Deep (60–100 m)'};

Temp0  = 0;
Month0 = 7;

nGrid = 200;

% File for saving/loading partial dependence results
pdFileName = fullfile( ...
    HGP.Data.Path2ResultsData, ...
    ['PD_ClustIDonTempDoyDepth' LM_Model '.mat']);

%% Calculate or load partial dependences

if calcPlot

    % Load trained model
    tmp = load(fullfile( ...
        HGP.Data.Path2ResultsData, ...
        ['trainedModel_ClustIDonTempDoyDepth' LM_Model '.mat']), ...
        'Mdl1');

    Mdl1 = tmp.Mdl1;
    clearvars tmp;

    % Define grid
    tempVec = linspace( ...
        min(tStPredictors.temp_imputed), ...
        max(tStPredictors.temp_imputed), ...
        nGrid)';

    monthVec = linspace( ...
        min(tStPredictors.Month_doy), ...
        max(tStPredictors.Month_doy), ...
        nGrid)';

    [tempGrid, monthGrid] = meshgrid(tempVec, monthVec);

    % Function returning class posterior probabilities
    
    fun = @(X) ml_returnPosterior(Mdl1, X);

    % Storage
    pds = cell(height(array2table(depths)), 1);
    xs  = cell(height(array2table(depths)), 1);
    ys  = cell(height(array2table(depths)), 1);

    for iD = 1:size(depths, 1)

        depthGrid = (depths(iD, 1):5:depths(iD, 2))';

        selectedData = table( ...
            repmat(Temp0,  size(depthGrid)), ...
            repmat(Month0, size(depthGrid)), ...
            depthGrid, ...
            'VariableNames', {'temp_imputed', 'Month_doy', 'Depth'});

        [pd, x, y] = partialDependence( ...
            fun, ...
            {'temp_imputed', 'Month_doy'}, ...
            selectedData, ...
            "QueryPoints", [tempGrid(:), monthGrid(:)], ...
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
        'Temp0', 'Month0', ...
        'tempVec', 'monthVec', ...
        'nGrid');

else

    % Load previously calculated partial dependences
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

    pcolor(x, y, classMap2D);

    clrs = f_Clrs_fresh(K);
    colormap(clrs)

    if iD == size(depths, 1)
        colorbar;
    end

    caxis([0.5, K + 0.5]);
    shading flat

    xlabel('Temperature');
    ylabel('Month');
    title(titles{iD});

    hold on;
    contour(x, y, pMax2D, ...
        [0, 0.4, 0.6, 0.8, 1], ...
        'k', ...
        'LineWidth', 1, ...
        'ShowText', 'on');

    %plot real data 
    indDepht = tStPredictors.Depth >= depths(iD,1)  & tStPredictors.Depth < depths(iD,2); 
    selectedPreds = tStPredictors(indDepht, :);
    scatter(selectedPreds.temp_imputed, selectedPreds.Month_doy, 20, tDMC.ClustID(indDepht), 'filled', 'MarkerEdgeColor', 'k');

    
 xData = selectedPreds.temp_imputed;
yData = selectedPreds.Month_doy;

% Remove NaNs
indGood = ~isnan(xData) & ~isnan(yData);
xData = xData(indGood);
yData = yData(indGood);

% Current plot limits, should correspond to the prediction grid
xLimits = [min(x(:)), max(x(:))];
yLimits = [min(y(:)), max(y(:))];

xlim(xLimits);
ylim(yLimits);

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

% Plot grey transparent mask
 plot(maskPoly, ...
     'FaceColor', [0.7 0.7 0.7], ...
     'FaceAlpha', 0.5, ...
     'EdgeColor', 'none');

    hold off;

end
AddLetters2Plots(gcf, {'d', 'e', 'f'}, HShift=-0.02, VShift=-0.1)