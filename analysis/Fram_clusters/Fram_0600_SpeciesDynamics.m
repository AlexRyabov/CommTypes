%% Plot the dynamics of species as a function of the day of year
proj_Initialize;

%define the grouping variables
day_step = 21; %step of the day of year 
tStPredictors.indDay = round(tStPredictors.doy/day_step)*day_step;
tStPredictors.indDepth = isdeep(tStPredictors.Depth);
tStPredictors.indEast =  tStPredictors.Lon > -0.01; %Grad0 station refers to east %  iseast(tStPredictors.Lon);

%find x most abundant families

%find mean abundance of each species across samples
% Find the histogram of average species abundance 
dSamples = Data.dSamplesNorm;  %species x samples
tSamples = Data.tSamplesNorm; %species x samples
idxSamplCols = DI.idxSamplCols;
tTraits = Data.tTraits;

%%find mean for each topGroupName in dSamples  (species x samples )
%1. Find sum for each group for each sample
topGroupName = 'Genus';
meanFamAbund = grpstats([tSamples(:, topGroupName), tSamples(:, idxSamplCols)],topGroupName, {'sum'});
meanFamAbund.GroupCount = [];
topGroupNames = meanFamAbund{:, topGroupName};
%transpose meanFamAbund to make statistics in samples
meanFamAbundSmpl = meanFamAbund{:, 2:end}';  %samples x topGroupName



%% 2. Find mean abundance  for given lon index, depth index and day index
%find mean abundance for each group
[meanDayDepthLonAbund, tDepthLonday] = grpstats(meanFamAbundSmpl, tStPredictors{:, {'indDepth', 'indEast', 'indDay'}}, {'mean', 'gname'});
%for control which samples we get 
meanDayDepthLonSamples = grpstats(tStPredictors(:, {'Lon', 'Lat', 'Depth', 'doy', 'indDepth', 'indEast', 'indDay'}), {'indDepth', 'indEast', 'indDay'}, {'mean', 'gname'});
%meanDayDepthLonAbund   depthlonday conmbinations x topGroupName (Genus)
tDepthLonday = cellfun(@str2double, tDepthLonday); %convert depth to double 
tDepthLonday =  array2table(tDepthLonday);
tDepthLonday.Properties.VariableNames = {'indDepth', 'indEast', 'indDay'};

meanDayDepthLonAbund = array2table(meanDayDepthLonAbund);
meanDayDepthLonAbund.Properties.VariableNames = topGroupNames;
% meanDayDepthLonAbund = movevars(meanDayDepthLonAbund, "Others");
topGroupNames = meanDayDepthLonAbund.Properties.VariableNames;
%save a temporal dataset



%% Filter shallow east

indCombinations = [ 0, 0; 0, 1; 1, 0; 1, 1];  %E/W; Shallow/Deep
Titles = [ "EGC, <30m", "EGC, \geq 30m", "WSC, <30", "WSC, \geq 30m"];
MaxTayaPerGroup =15;
AllSpUnique = {};

% Calculate the thirds of observed period
% dayRange = range(TDayFiltered.indDay);
% firstThird = min(TDayFiltered.indDay) + dayRange / 3;
% secondThird = min(TDayFiltered.indDay) + 2 * dayRange / 3;
DayRanges = linspace(min(tDepthLonday.indDay), max(tDepthLonday.indDay), 4);  %define the boundaries 
for iInds = 1:height(indCombinations )
    %select date for given location and depth 
    ind = tDepthLonday.indEast == indCombinations(iInds, 1) &  tDepthLonday.indDepth == indCombinations(iInds, 2) ;
    %get groups of spring, summer, and autumn species
    [groupedTables{iInds}, AllSpecies] = findSuccessionByGroups(meanDayDepthLonAbund(ind, :) ,  tDepthLonday(ind, :), DayRanges,  MaxTayaPerGroup);
    AllSpUnique = [AllSpUnique, AllSpecies];
end
AllSpUnique = unique(AllSpUnique)';  %create a list of unique species identifications (genera, class, phyla etc)


        %assign colors automatically
        [DictsColrs] = colorDict(AllSpUnique);

%%
fg1 = f_MakeFigure(1, [-1, -1, 1.4613e+03 741]);
clf
%create tiled layout
%tiledlayout(3, 4, "TileSpacing","tight", "Padding","tight");
tiledlayout(3, 4, "TileSpacing","compact", "Padding","tight");
PlotAxes = arrayfun(@(x) nexttile, 1:12, 'UniformOutput', false);
for iInds = 1:height(indCombinations )
    %define index of subplots
    ind = tDepthLonday.indDepth == indCombinations(iInds, 2) & tDepthLonday.indEast == indCombinations(iInds, 1);
    PlotInds = [iInds, iInds+4, iInds+8];
    %plot results 
    AddYLabel = iInds==1; 
    plotSpeciesGroups(Titles(iInds), PlotAxes(PlotInds), groupedTables{iInds}, tDepthLonday(ind, :), DictsColrs, AddYLabel)
end
AddLetters2Plots('Location',   'SouthWest', 'VShift', -0.03)


