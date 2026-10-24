% Plot environmental conditions
proj_Initialize;

%% Imputation statistics

imputationFields = ["temp"; "DIN"; "DIN_100_200"; ...
                    "nitrate"; "phosphate"; "silicate"];
imputedFields = imputationFields + "_imputed";
nSamples = height(tStPredictors);

imputationStats = table(imputationFields, repmat(nSamples, numel(imputationFields), 1), ...
                        zeros(numel(imputationFields), 1), ...
                        zeros(numel(imputationFields), 1), ...
                        zeros(numel(imputationFields), 1), ...
                        zeros(numel(imputationFields), 1), ...
                        'VariableNames', {'Field', 'NSamples', 'OriginalMissingN', ...
                                          'ImputedColumnNonMissingN', 'FilledMissingN', ...
                                          'ImputedPercent'});

for iF = 1:numel(imputationFields)
    originalMissing = ismissing(tStPredictors.(imputationFields(iF)));
    imputedPresent = ~ismissing(tStPredictors.(imputedFields(iF)));

    imputationStats.OriginalMissingN(iF) = sum(originalMissing);
    imputationStats.ImputedColumnNonMissingN(iF) = sum(imputedPresent);
    imputationStats.FilledMissingN(iF) = sum(originalMissing & imputedPresent);
    imputationStats.ImputedPercent(iF) = imputationStats.FilledMissingN(iF) / nSamples * 100;
end

disp('Imputation statistics for environmental predictors:')
disp(imputationStats)

[fg, MaxFigure] = addFindFigure('envCond', MaxFigure, ...
                                [-1, -1, 1222 351], 'normal', 'on');
clf

tl = tiledlayout(1, 3);
tl.TileSpacing = "compact";
tl.Padding = "compact";


%% a: Seasonal decline in surface PAR

ax1 = nexttile(tl, 1);

nYears = max(tStPredictors.Year) - min(tStPredictors.Year) + 1;
colormap(ax1, jet(nYears));

scatter(ax1, tStPredictors.Month_doy, tStPredictors.mod_aq_par, ...
        30, tStPredictors.Year, 'filled');

clim(ax1, [min(tStPredictors.Year)-0.5, max(tStPredictors.Year)+0.5]);
xlim(ax1, [5.9, 10.1]);

xticks(ax1, 6:10);
xticklabels(ax1, {'Jun', 'Jul', 'Aug', 'Sep', 'Oct'});

xlabel(ax1, 'Month');
ylabel(ax1, 'Surface PAR (\mumol photons m^{-2} s^{-1})');

lm = stepwiselm(tStPredictors.Month_doy, tStPredictors.mod_aq_par, 'quadratic');
xVals = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 100);
yVals = predict(lm, xVals');

hold(ax1, 'on');
plot(ax1, xVals, yVals, 'k', 'LineWidth', 1.2);
colormap(ax1, "jet");

cb1 = colorbar(ax1);
cb1.Label.String = 'Year';

box(ax1, 'on');


%% b: Nitrate depth profile

ax2 = nexttile(tl, 2);

yJitter = (rand(size(tStPredictors.Depth)) - 0.5) * 2;

s = scatter(ax2, tStPredictors.nitrate, tStPredictors.Depth + yJitter, ...
            25, 'filled');
s.MarkerFaceAlpha = 0.3;

set(ax2, 'YDir', 'reverse');
set(ax2, 'XScale', 'log');

ylim(ax2, [0, 100]);
xlim(ax2, [1e-2, 20]);

xlabel(ax2, 'Nitrate (\mumol L^{-1})');
ylabel(ax2, 'Depth (m)');

box(ax2, 'on');


%% c: Temperature across longitude and depth

ax3 = nexttile(tl, 3);

s = scatter(ax3, tStPredictors.Lon, tStPredictors.Depth + yJitter, ...
            30, tStPredictors.temp, 'filled');
s.MarkerFaceAlpha = 0.5;

colormap(ax3, "jet");

set(ax3, 'YDir', 'reverse');
ylim(ax3, [0, 100]);

xlabel(ax3, 'Longitude');
ylabel(ax3, 'Depth (m)');

cb3 = colorbar(ax3);
cb3.Label.String = 'Temperature (^{\circ}C)';

box(ax3, 'on');


%% Panel letters

AddLetters2Plots(fg, {'a', 'b', 'c'}, ...
                 'HShift', -0.05, ...
                 'Location', 'NorthWest');
