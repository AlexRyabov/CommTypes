%plot environmental conditions 
proj_Initialize;

%%
[fg, MaxFigure] = addFindFigure('envCond', MaxFigure, [-1, -1,  824.6667 604], 'normal', 'on');
t = nexttile
nYears = max(tStPredictors.Year) - min(tStPredictors.Year) + 1;
colormap(t, distinguishable_colors(nYears))
scatter(tStPredictors.Month_doy, tStPredictors.mod_aq_par, 30, tStPredictors.Year, 'filled');
clim([min(tStPredictors.Year)-0.5, max(tStPredictors.Year) + 0.5])
xlim([5.9, 10.1])
xlabel('Month')
ylabel('PAR,  µmol photons m^{-2} s^{-1}');
lm = stepwiselm(tStPredictors.Month_doy, tStPredictors.mod_aq_par, 'quadratic');
xVals = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 20);
yVals = lm.predict(xVals');
hold on
plot(xVals, yVals, 'k')
cb = colorbar
cb.Label.String = 'Year';

t = nexttile
scatter(tStPredictors.nitrate, tStPredictors.Depth, 'filled');
xscale('log')
set(gca,'YDir','reverse')
ylim([0, 50])
xlabel('Nitrate')
ylabel('depth, m')
% ylabel('PAR,  µmol photons m^{-2} s^{-1}');
% lm = stepwiselm(tStPredictors.Month_doy, tStPredictors.mod_aq_par, 'quadratic');
% xVals = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 20);
% yVals = lm.predict(xVals');
% hold on
% plot(xVals, yVals, 'k')


t = nexttile
scatter(tStPredictors.Lon, tStPredictors.Depth, 30, tStPredictors.temp,   'filled');
colormap(t, "jet")
set(gca,'YDir','reverse')
cb = colorbar
cb.Label.String = 'T, °C';    
ylim([0, 100])
xlabel('Longitude')
ylabel('depth, m')
% ylabel('PAR,  µmol photons m^{-2} s^{-1}');
% lm = stepwiselm(tStPredictors.Month_doy, tStPredictors.mod_aq_par, 'quadratic');
% xVals = linspace(min(tStPredictors.Month_doy), max(tStPredictors.Month_doy), 20);
% yVals = lm.predict(xVals');
% hold on
% plot(xVals, yVals, 'k')
AddLetters2Plots('HShift', 0.0,  'Location', 'SouthWest', 'FitLocation', true)
