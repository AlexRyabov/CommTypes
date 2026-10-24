function [ax, SampleStat] = violLogWithNames(Samples, SamplesMean, face_color, ylims, varargin)
%
p = inputParser;
p.addParameter('Rotate',45,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('YScale','log',@(x)ischar(x)||isstring(x));
p.addParameter('YLabel','',@(x)ischar(x)||isstring(x));
p.addParameter('YQuantiles',[],@(x)isnumeric(x));
p.parse(varargin{:});
tickRotation = p.Results.Rotate;
yl  = string(p.Results.YLabel);
threshold = log10(1.05);

if size(face_color, 1) == 1
    face_colors = repmat(face_color, height(Samples), 1);
elseif size(face_color, 1) == height(Samples) && size(face_color, 2) == 3
    face_colors = face_color;
else
    error('face_color must be a 1-by-3 RGB vector or a height(Samples)-by-3 RGB matrix.');
end

%% make a paler version of your base color (toward white)

lighten = @(c,a) c + (1-c).*a;     % a in [0,1], larger = paler
SampleStat.val_std = NaN(height(Samples), 1);
SampleStat.log_norm_val_mean = SampleStat.val_std;
SampleStat.log_norm_val_median = SampleStat.val_std;
for iR = 1:height(Samples)
    vals = Samples{Samples.Row(iR), :};
    pos = vals > 0 & ~isnan(vals);
    if strcmp(p.Results.YScale, 'log')
        if any(~pos)
            epsVal = max(1e-6, min(vals(pos))/10);
            vals(~pos) = epsVal;
        end
        vals = log10(vals/SamplesMean{Samples.Row(iR), "mean_diversity"});
    else
        vals = (vals/SamplesMean{Samples.Row(iR), "mean_diversity"});
    end

    SampleStat.val_std(iR) = std(vals);
    SampleStat.log_norm_val_mean(iR) = mean(vals);
    SampleStat.log_norm_val_median(iR) = median(vals);
    row_color = face_colors(iR, :);
    pale_color = lighten(row_color, 0.6);
    if SampleStat.log_norm_val_median(iR) >  threshold
        viol_alpha = 0.5;
        clr = row_color;
    else
        viol_alpha = 0.25;
        clr = pale_color;
    end

    vs = Violin({vals},iR, 'ViolinColor', {clr}, 'ViolinAlpha', {viol_alpha}, 'ShowData', false, 'MedianColor', row_color, ...
        'QuartileStyle','shadow', 'MedianMarkerSize', 12);

end
%%

x = 1:1:height(Samples);
 small_val_ind = SampleStat.log_norm_val_median <  threshold;
ax = gca;
ylim(ax, ylims)
ax.XTick = x;
ax.FontSize = 8;

names = Samples.Row;
SampleStat.names = names;
ax = addTickLabels(names,ax,small_val_ind,tickRotation);

ylabel(ax, yl);
hold on
yline( 0, 'k');
hold off
end
