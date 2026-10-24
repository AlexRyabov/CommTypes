function fgNum = hg_cd_aux_plotspeceis(SpeciesAbund, fgName, fgNum, color)
arguments
    SpeciesAbund   %samples (rows) x species (columns)
    fgName
    fgNum
    color = [];
end

%%
meanAbund = mean(SpeciesAbund, 1);


[~, ind] = sort(meanAbund{1, :}, 'descend');
[fg, fgNum] = addFindFigure(fgName, fgNum, [], 'docked', 'on');

%select first 100 elements to plot

if isempty(color)
    m = min(100, width(meanAbund));
    bar(meanAbund{1, ind(1:m)})
    ax = gca; ax.YScale = 'log';
else
    m = width(meanAbund);
    x = 1:m;
    meanABund2plot = meanAbund{1, ind(1:m)};
    t = nexttile(2)
    scatter(x, meanABund2plot,(color(ind(1:m))+0.01).^2*100, color(ind(1:m)), 'filled');
    t.YScale = 'log';
    cb = colorbar;
    cb.Label.String = 'Importance';
    xlabel('rank'); ylabel('Abundance');
    t = nexttile(1)
    [mostImp50, mostImpInd]= maxk(color(ind(1:m)), 50);
    scatter(color(ind(1:m)), meanABund2plot, 'filled');
    hold on
    scatter(mostImp50, meanABund2plot(mostImpInd), 'filled');
    hold off
    t.YScale = 'log';
    xlabel('Importance'); 
end    
ax = gca; ax.YScale = 'log';
end