
 data = load(fullfile(HGP.Data.Path2ResultsData, ['trainedLocationModel_' LM_Model '.mat']));
% 
 Mdl1 = data.Mdl1;
 ClassifID = data.ClassifID;
 mask_seq2 = data.mask_seq2;
 cv_inner = data.cv_inner;
 bestParams = data.bestParams;
 cvErr = data.cvErr;
 ml_train_func = @ml_train_ECOC;
% 

%% INPUTS
% Mdl1   : trained TreeBagger classifier (Method='classification')
% X_data : table with ONLY numeric predictor columns used to train Mdl1
%         (same variable names as Mdl1.PredictorNames)
%rng(42);

predNames  = Mdl1.PredictorNames;




if isprop(Mdl1,'OOBPermutedPredictorDeltaError')
    pred_importance = Mdl1.OOBPermutedPredictorDeltaError;
else
   [pred_importance, impTable, baseLoss] = permutationImportance(X_train(:, predNames), Y_train, Mdl1, cv_inner);
    % [~, idx] = sort(pred_importance, "descend");    
    %tImp =  permutationImportance(Mdl1);
    %pred_importance = tImp.ImportanceMean;
end
[~, idx] = sort(pred_importance, "descend");
predNames = predNames(idx);

[fg, MaxFigure] = addFindFigure('Importance', MaxFigure, ...
    [10, 20, 400, 400], 'normal', 'on');
bar( predNames, pred_importance(idx))




%%
calcPD = true;    % true = compute and save, false = load and plot only
classNames = string(Mdl1.ClassNames);

K      = numel(classNames);      % number of classes
nPred  = numel(predNames);       % number of predictors
nGrid  = 20;                     % number of grid points on xs
B      = 0;                     % number of bootstrap replicates

Xpd = X_train(:, predNames);     % table with predictors

fun = @(X) ml_returnPosterior(Mdl1, X);

% File to store partial dependence results
pdFile = fullfile(HGP.Data.Path2ResultsData, ['PDP_results_' LM_Model '.mat']) ;

% Storage for PDP and CIs for each predictor
pdMeanAll  = cell(nPred, 1);   % each: [K x nGrid]
pdLowAll   = cell(nPred, 1);   % each: [K x nGrid]
pdHighAll  = cell(nPred, 1);   % each: [K x nGrid]
xsAll      = cell(nPred, 1);   % each: [1 x nGrid]

N = height(Xpd);

if calcPD
    % =========================
    % First loop: compute PDPs + CIs
    % =========================
    for j = 1:nPred
        predj = predNames{j}

        % Build grid for this predictor
        xj   = Xpd.(predj);
        xmin = min(xj);
        xmax = max(xj);
        xs   = linspace(xmin, xmax, nGrid);   % 1 x nGrid

        % 1) PDP on full data (mean curve)
        pdFull = partialDependence(fun, predj, Xpd, ...
            "UseParallel",   true, ...
            "QueryPoints",   xs', ...          % xs': nGrid x 1
            "OutputColumns", 1:K);             % pdFull: K x nGrid

        % 2) Bootstrap distribution of PDPs
        pdBoot = NaN(B, K, nGrid);             % [B x K x nGrid]

        parfor b = 1:B
            idx = randi(N, N, 1);              % indices with replacement
            Xb  = Xpd(idx, :);                 % bootstrap sample of data
            Yb     = Y_train(idx);
            switch ClassifID
                case 'ECOC_SVM'
                    Mdl_b = ml_train_func(Xb, Yb, bestParams,  'FitPosterior', true);
                otherwise
                    Mdl_b = ml_train_func(Xb, Yb, bestParams);
            end

            funB = @(X) ml_returnPosterior(Mdl_b, X);
            pd_b = partialDependence(funB, predj, Xb, ...
                "UseParallel",   false, ...
                "QueryPoints",   xs', ...
                "OutputColumns", 1:K);         % [K x nGrid]

            pdBoot(b, :, :) = pd_b;
        end

        % 3) Compute pointwise CIs over bootstrap dimension (1st dim)
        %    Result: K x nGrid
        pdLow  = squeeze(prctile(pdBoot,  2.5, 1));
        pdHigh = squeeze(prctile(pdBoot, 97.5, 1));

        % Store
        pdMeanAll{j} = pdFull;
        pdLowAll{j}  = pdLow;
        pdHighAll{j} = pdHigh;
        xsAll{j}     = xs;      % keep row vector [1 x nGrid]
    end

    % Pack everything into a struct for saving
    PD = struct();
    PD.pdMeanAll   = pdMeanAll;
    PD.pdLowAll    = pdLowAll;
    PD.pdHighAll   = pdHighAll;
    PD.xsAll       = xsAll;
    PD.predNames   = predNames;
    PD.classNames  = classNames;
    PD.nGrid       = nGrid;
    PD.B           = B;
    PD.K           = K;

    save(pdFile, 'PD');   % save all PDP results

else
    % =========================
    % Load precomputed PDPs
    % =========================
    S = load(pdFile, 'PD');
    PD = S.PD;

    % Restore variables from struct (optional, для удобства)
    pdMeanAll  = PD.pdMeanAll;
    pdLowAll   = PD.pdLowAll;
    pdHighAll  = PD.pdHighAll;
    xsAll      = PD.xsAll;
    predNames  = PD.predNames;
    classNames = PD.classNames;
    nGrid      = PD.nGrid;
    B          = PD.B;
    K          = PD.K;
end

%% =========================
%  Second loop: plotting
%  ==========================
[fg, MaxFigure] = addFindFigure('Part depend', MaxFigure, ...
    [10, 20, 1145, 656], 'normal', 'on');

tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
clrs = f_Clrs_fresh(K);

for j = 1:nPred
    predj = predNames{j};
    pd    = pdMeanAll{j};      % [K x nGrid]
    pdLow = pdLowAll{j};       % [K x nGrid]
    pdHigh= pdHighAll{j};      % [K x nGrid]
    xs    = xsAll{j};          % [1 x nGrid]

    % Range of mean PDP per class
    pdRange = max(pd, [], 2) - min(pd, [], 2);   % K x 1

    % Map range in [0,1] to line width in [1,10]
    lineWidths = 1 + 9 * pdRange;                % K x 1

    ax = nexttile;
    hold(ax, 'on');
    colororder(ax, clrs);

    for k = 1:K
        % Confidence band for class k
        xPatch = [xs, fliplr(xs)];
        yPatch = [pdLow(k, :), fliplr(pdHigh(k, :))];

        p = patch('XData', xPatch, 'YData', yPatch, ...
              'FaceColor', clrs(k, :), ...
              'FaceAlpha', 0.15, ...
              'EdgeColor', 'none', ...
              'Parent', ax);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        % Mean curve
        plot(ax, xs, pd(k, :), ...
            'LineWidth', lineWidths(k), ...
            'Color', clrs(k, :));
    end

    xlabel(ax, predj, 'Interpreter','none');
    ylabel(ax, 'Partial dependence: P(class | x)');
    title(ax, sprintf('PDP: %s', predj), 'Interpreter','none');

    hold(ax, 'off');
    drawnow;
end

% One legend for all tiles
lg = legend(classNames, 'NumColumns', min(3, K), 'Location', 'southoutside');
lg.Title.String = 'Classes';

