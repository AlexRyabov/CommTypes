function ax = barLogWithNames(vals, names, face_color, ylims, varargin)
% vals: numeric vector
% names: string/cellstr vector, same length as vals

p = inputParser;
p.addParameter('Rotate',45,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('YLabel','',@(x)ischar(x)||isstring(x));
p.addParameter('YQuantiles',[],@(x)isnumeric(x));
p.parse(varargin{:});
rot = p.Results.Rotate;
yl  = string(p.Results.YLabel);

pos = vals > 0 & ~isnan(vals);
if any(~pos)
    epsVal = max(1e-6, min(vals(pos))/10);
    vals(~pos) = epsVal;
end
ax = gca;
%b = bar(ax, vals, 'FaceColor',face_color);
x = 1:length(vals);
b = bar(x, vals, 'FaceColor','flat', 'EdgeColor','none');
if ~isempty(p.Results.YQuantiles)
    hold on
    errorbar(x, vals, vals-p.Results.YQuantiles(:, 1), vals+p.Results.YQuantiles(:, 2), 'Color', 'k')
    hold off
end
% make a paler version of your base color (toward white)
lighten = @(c,a) c + (1-c).*a;     % a in [0,1], larger = paler
pale_color = lighten(face_color, 0.6);
% assign colors per bar
C = repmat(face_color, numel(vals), 1);
small_val_ind = vals < 1.05;
C(small_val_ind, :) = repmat(pale_color, nnz(small_val_ind), 1);
b.CData = C;

% set(ax, 'YScale','log'); %grid(ax,'on'); ax.YMinorGrid='on';
% ylim(ylims)
% ax.XTick = 1:numel(vals);
% if ~isempty(names)
%     ax.XTickLabel = string(names);
%     ax.XTickLabelRotation = rot;
%     ax.TickLabelInterpreter = 'none';
% else
%     ax.XTickLabel = string(names);
% end
ax.YScale = 'log';
ylim(ax, ylims)
ax.XTick = 1:numel(vals);
ax.FontSize = 8;
if ~isempty(names)
    labs = string(names);

    % Use TeX markup and escape underscores
    ax.TickLabelInterpreter = 'tex';
    labs = replace(labs, "_", "\_");

    % Bold only where vals > 1
    mask = ~small_val_ind;
    labs(mask) = "{\bf " + labs(mask) + "}";

    fsSmall  = max(6, round(ax.FontSize*0.75));  % small font
    smallMask = ~mask;
    labs(smallMask) = "{\fontsize{" + string(fsSmall) + "} " + labs(smallMask) + "}";

    ax.XTickLabel = labs;
    ax.XTickLabelRotation = rot;
else
    ax.XTickLabel = string(names);
end

ylabel(ax, yl);
hold on
plot( [1, numel(vals)], [1, 1], 'k');
hold off
end
