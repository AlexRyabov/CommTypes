function ax = addTickLabels(names,ax,small_val_ind,rot)
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
end