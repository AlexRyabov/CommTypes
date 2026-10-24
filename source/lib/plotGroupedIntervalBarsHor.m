function [fg, ax, cb, maxAbsValue] = plotGroupedIntervalBarsHor(values, valuesLow, valuesHigh, xTickLabels, yTickLabels, yGroups, yRoles, opts)
% plotGroupedIntervalBarsHor Plot horizontal signed bars with intervals.
%
% values, valuesLow, and valuesHigh are [nY x nX] matrices. Bars are drawn
% horizontally inside each x group, with yGroups shown as right-side labels
% and separators. yRoles can be used to color y tick labels and draw a role
% RoleColors color y tick labels. No separate role legend is drawn.

    if nargin < 8 || isempty(opts)
        opts = struct();
    end

    values = double(values);
    valuesLow = double(valuesLow);
    valuesHigh = double(valuesHigh);
    xTickLabels = string(xTickLabels(:))';
    yTickLabels = string(yTickLabels(:));
    yGroups = string(yGroups(:));
    yRoles = string(yRoles(:));

    [nY, nX] = size(values);

    if ~isequal(size(valuesLow), [nY, nX]) || ~isequal(size(valuesHigh), [nY, nX])
        error('values, valuesLow, and valuesHigh must have the same size.');
    end

    if numel(xTickLabels) ~= nX
        error('xTickLabels must have one label per column in values.');
    end

    if numel(yTickLabels) ~= nY
        error('yTickLabels must have one label per row in values.');
    end

    if isempty(yGroups)
        yGroups = strings(nY, 1);
    elseif numel(yGroups) ~= nY
        error('yGroups must be empty or have one label per row in values.');
    end

    if isempty(yRoles)
        yRoles = repmat("Unknown", nY, 1);
    elseif numel(yRoles) ~= nY
        error('yRoles must be empty or have one label per row in values.');
    end

    fg = local_getOpt(opts, 'Figure', []);
    if isempty(fg) || ~isgraphics(fg, 'figure')
        figureName = local_getOpt(opts, 'FigureName', 'Grouped horizontal interval bars');
        figurePosition = local_getOpt(opts, 'FigurePosition', [24 332 1461 455]);
        fg = figure('Name', figureName, 'Position', figurePosition, 'Color', 'w');
    else
        fg.Color = 'w';
    end

    if local_getOpt(opts, 'ClearFigure', true)
        clf(fg);
    end

    ax = local_getOpt(opts, 'Axes', []);
    if isempty(ax) || ~isgraphics(ax, 'axes')
        ax = axes(fg);
    end

    hold(ax, 'on');

    barHalfWidth = local_getOpt(opts, 'BarHalfWidth', 0.42);
    barHalfHeight = local_getOpt(opts, 'BarHalfHeight', 0.32);
    intervalCapHalfHeight = local_getOpt(opts, 'IntervalCapHalfHeight', 0.10);
    cmap = local_getOpt(opts, 'Colormap', local_blueWhiteRed(256));

    maxAbsValue = local_getOpt(opts, 'MaxAbsValue', []);

    if isempty(maxAbsValue)
        plotValues = abs([values(:); valuesLow(:); valuesHigh(:)]);
        plotValues = plotValues(plotValues > 0 & isfinite(plotValues));

        if isempty(plotValues)
            maxAbsValue = 1;
        else
            maxAbsValue = max(plotValues);
        end
    end

    if maxAbsValue == 0 || ~isfinite(maxAbsValue)
        maxAbsValue = 1;
    end

    for iY = 1:nY
        for iX = 1:nX
            value = values(iY, iX);
            valueLow = valuesLow(iY, iX);
            valueHigh = valuesHigh(iY, iX);

            if ~isfinite(value) || ~isfinite(valueLow) || ~isfinite(valueHigh)
                continue
            end

            x0 = iX;
            xMean = iX + barHalfWidth * value / maxAbsValue;
            xLow = iX + barHalfWidth * valueLow / maxAbsValue;
            xHigh = iX + barHalfWidth * valueHigh / maxAbsValue;
            xLo = min(xLow, xHigh);
            xHi = max(xLow, xHigh);
            clr = local_valueColor(value, maxAbsValue, cmap);

            patch(ax, ...
                [x0 xMean xMean x0], ...
                [iY - barHalfHeight, iY - barHalfHeight, ...
                 iY + barHalfHeight, iY + barHalfHeight], ...
                clr, ...
                'EdgeColor', 'none', ...
                'FaceAlpha', local_getOpt(opts, 'FaceAlpha', 0.92));
            plot(ax, [xLo xHi], [iY iY], '-', ...
                'Color', local_getOpt(opts, 'IntervalColor', [0 0 0]), ...
                'LineWidth', local_getOpt(opts, 'IntervalLineWidth', 1.05));
            plot(ax, [xLo xLo], [iY - intervalCapHalfHeight, iY + intervalCapHalfHeight], '-', ...
                'Color', local_getOpt(opts, 'IntervalColor', [0 0 0]), ...
                'LineWidth', local_getOpt(opts, 'IntervalCapLineWidth', 0.85));
            plot(ax, [xHi xHi], [iY - intervalCapHalfHeight, iY + intervalCapHalfHeight], '-', ...
                'Color', local_getOpt(opts, 'IntervalColor', [0 0 0]), ...
                'LineWidth', local_getOpt(opts, 'IntervalCapLineWidth', 0.85));
        end
    end

    for iX = 1:nX
        xline(ax, iX, '-', ...
            'Color', local_getOpt(opts, 'XGridLineColor', [0.75 0.75 0.75]), ...
            'LineWidth', local_getOpt(opts, 'XGridLineWidth', 0.5));
    end

    set(ax, 'YDir', 'reverse');
    xlim(ax, [0.5 nX + 1.25]);
    ylim(ax, [0.5 nY + 0.5]);

    xticks(ax, 1:nX);
    xticklabels(ax, xTickLabels);
    yticks(ax, 1:nY);

    roleColors = local_getOpt(opts, 'RoleColors', []);

    if isempty(roleColors)
        yticklabels(ax, cellstr(yTickLabels));
        ax.TickLabelInterpreter = 'none';
    else
        yticklabels(ax, local_coloredYTickLabels(yTickLabels, yRoles, roleColors));
        ax.TickLabelInterpreter = 'tex';
    end

    ax.FontSize = local_getOpt(opts, 'FontSize', 8);
    xlabel(ax, local_getOpt(opts, 'XLabel', ''));
    yLabelHandle = ylabel(ax, local_getOpt(opts, 'YLabel', ''));
    title(ax, local_getOpt(opts, 'Title', ''));

    yLabelHandle.Units = 'normalized';
    yLabelHandle.Position(1) = -0.19;
    yLabelHandle.Position(2) = 0.5;

    grid(ax, 'off');
    box(ax, 'on');
    colormap(ax, cmap);
    clim(ax, [-maxAbsValue maxAbsValue]);
    cb = colorbar(ax);
    cb.Label.String = local_getOpt(opts, 'ColorbarLabel', '');

    if local_getOpt(opts, 'ShowGroupLabels', true) && any(strlength(yGroups) > 0)
        local_addYGroupLabels(ax, yGroups, nX, local_getOpt(opts, 'GroupFontSize', ax.FontSize));
    end

end

function value = local_getOpt(opts, name, defaultValue)

    if isstruct(opts) && isfield(opts, name)
        value = opts.(name);
    else
        value = defaultValue;
    end
end

function local_addYGroupLabels(ax, yGroups, nX, fontSize)

    for ii = 1:(numel(yGroups) - 1)
        if yGroups(ii) ~= yGroups(ii + 1)
            yline(ax, ii + 0.5, '-', ...
                'Color', [0.25 0.25 0.25], ...
                'LineWidth', 0.9);
        end
    end

    uniqueGroups = unique(yGroups, 'stable');

    for iGroup = 1:numel(uniqueGroups)
        idxGroup = find(yGroups == uniqueGroups(iGroup));
        yMid = mean(idxGroup);

        text(ax, nX + 0.23, yMid, uniqueGroups(iGroup), ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', fontSize, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end
end

function tickLabels = local_coloredYTickLabels(labels, roleNames, roleColors)

    labels = string(labels);
    roleNames = string(roleNames);
    tickLabels = cell(numel(labels), 1);

    for i = 1:numel(labels)
        clr = local_roleColor(roleNames(i), roleColors);
        label = local_escapeTex(labels(i));

        tickLabels{i} = sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%s', ...
            clr(1), clr(2), clr(3), label);
    end
end

function clr = local_roleColor(roleName, roleColors)

    roleName = char(string(roleName));

    if isa(roleColors, 'containers.Map') && isKey(roleColors, roleName)
        clr = roleColors(roleName);
    elseif isa(roleColors, 'containers.Map') && isKey(roleColors, 'Unknown')
        clr = roleColors('Unknown');
    else
        clr = [0.15 0.15 0.15];
    end
end

function label = local_escapeTex(label)

    label = char(string(label));
    label = strrep(label, '\', '\\');
    label = strrep(label, '_', '\_');
    label = strrep(label, '^', '\^');
    label = strrep(label, '{', '\{');
    label = strrep(label, '}', '\}');
end

function cmap = local_blueWhiteRed(n)

    if nargin < 1
        n = 256;
    end

    n1 = floor(n / 2);
    n2 = n - n1;

    blue = [0.20 0.38 0.80];
    white = [1.00 1.00 1.00];
    red = [0.82 0.18 0.12];

    cmap1 = [linspace(blue(1), white(1), n1)', ...
             linspace(blue(2), white(2), n1)', ...
             linspace(blue(3), white(3), n1)'];

    cmap2 = [linspace(white(1), red(1), n2)', ...
             linspace(white(2), red(2), n2)', ...
             linspace(white(3), red(3), n2)'];

    cmap = [cmap1; cmap2];
end

function clr = local_valueColor(value, maxAbsValue, cmap)

    if maxAbsValue <= 0 || ~isfinite(maxAbsValue)
        maxAbsValue = 1;
    end

    t = (value / maxAbsValue + 1) / 2;
    t = max(0, min(1, t));
    idx = 1 + round(t * (size(cmap, 1) - 1));
    clr = cmap(idx, :);
end
