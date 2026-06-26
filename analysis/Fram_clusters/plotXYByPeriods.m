function plotXYByPeriods(edges, names, split_value, XY, groupID,  varargin)
% plotXYByPeriods
%   edges     : vector of bin edges for split_value (e.g., [6 6.5 7 7.5 8 9 10])
%   names     : cellstr of panel titles, length(names) = numel(edges)-1
%   XY        : Nx2 numeric [X=Lon, Y=Depth]
%   groupID   : Nx1 labels (numeric/string/categorical/cellstr)
%   split_value : Nx1 numeric (e.g., 6.0..9.99)
%
% Optional name-value pairs:
%   'MarkerSize' (40), 'LonJitter' (0.2), 'DepthJitter' (3), 'Alpha' (0.7)
%   'XLabel' ('Longitude'), 'YLabel' ('Depth (m)')

    % ---- options
    p = inputParser;
    p.addParameter('MarkerSize', 40, @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('LonJitter',  0.2, @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('DepthJitter', 3,  @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('Alpha',     1,  @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('XLabel', 'Longitude', @(x)ischar(x)||isstring(x));
    p.addParameter('YLabel', 'Depth (m)', @(x)ischar(x)||isstring(x));
    p.parse(varargin{:});
    mark_size    = p.Results.MarkerSize;
    lon_jitter   = p.Results.LonJitter;
    depth_jitter = p.Results.DepthJitter;
    alphaVal     = p.Results.Alpha;
    xlab         = string(p.Results.XLabel);
    ylab         = string(p.Results.YLabel);

    % ---- checks
    if size(XY,2) ~= 2
        error('XY must be an Nx2 matrix: [:,1]=X (Lon), [:,2]=Y (Depth).');
    end
    if numel(names) ~= numel(edges)-1
        error('Length of names must be numel(edges)-1.');
    end
    N = size(XY,1);
    if ~isequal(N, numel(groupID), numel(split_value))
        error('XY, groupID, and split_value must have the same number of rows/elements.');
    end

    % ---- groups & colors
    grpStr = string(groupID);                      % robust to numeric/cellstr/categorical
    labels = unique(grpStr, 'sorted');             % group labels as strings
    K      = numel(labels);
    clrs = f_Clrs_fresh(K);

    % ---- layout
    nP   = numel(names);
 
  
    tL   = tiledlayout('TileSpacing','compact','Padding','compact');

    % global limits for consistency across panels
    X = XY(:,1); Y = XY(:,2); M = split_value(:);
    goodAll = ~isnan(X) & ~isnan(Y) & ~isnan(M);
    xlimAll = [min(X(goodAll),[],'omitnan'), max(X(goodAll),[],'omitnan')];
    ylimAll = [min(Y(goodAll),[],'omitnan'), max(Y(goodAll),[],'omitnan')];

    for i = 1:nP
        ax = nexttile(tL, i); hold(ax,'on');

        % Bin: [edges(i), edges(i+1)) except last bin inclusive on right
        if i < nP
            inBin = (M >= edges(i)) & (M < edges(i+1));
        else
            inBin = (M >= edges(i)) & (M <= edges(i+1));
        end

        % one handle per group for legend
        h = gobjects(K,1);

        for c = 1:K
            inGrp = (grpStr == labels(c));
            idx   = inBin & inGrp & goodAll;

            n = sum(idx);
            if n > 0
                jX = (rand(n,1)*2 - 1) * lon_jitter;
                jY = (rand(n,1)*2 - 1) * depth_jitter;

                h(c) = scatter(ax, X(idx) + jX, Y(idx) + jY, ...
                    mark_size, repmat(clrs(c,:), n, 1), ...
                    'filled', 'MarkerEdgeColor','k', ...
                    'DisplayName', sprintf('C %s', labels(c)));

                h(c).MarkerFaceAlpha = alphaVal;
            end
        end

        set(ax,'YDir','reverse');    % depth increases downward
        grid(ax,'on');
        xlabel(ax, xlab);
        ylabel(ax, ylab);
        title(ax, names{i});

        xlim(ax, xlimAll);
        ylim(ax, ylimAll);

        valid = h(isgraphics(h));
        if ~isempty(valid)
            legend(ax, valid, 'Location','bestoutside');
        end
    end
end
