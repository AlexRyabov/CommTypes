function out = gmm_cluster(tDMC, varargin)
%GMM_CLUSTER_TDMC  Gaussian mixture clustering for Diffusion Map coordinates.
% INPUT:
%   tDMC : table, N×q with numeric DM coordinates (no normalization is applied).
%
% NAME-VALUE PAIRS (optional):
%   'KRange'           : vector of candidate components, default 1:10
%   'CovarianceTypes'  : cellstr of {'full','diag'}, default {'full','diag'}
%   'SharedCovariance' : logical vector [false true] to try both, default [false true]
%   'Replicates'       : integer, fitgmdist replicates, default 10
%   'Regularization'   : scalar regularization value, default 1e-6
%   'MaxIter'          : max EM iterations, default 1000
%   'PlotDims'         : 1×2 indices of columns to plot (if available), default [1 2]
%   'RandomState'      : RNG seed (scalar) or [] to skip, default 42
%
% OUTPUT (struct):
%   out.gm              : best fitted gmdistribution
%   out.P               : N×K posterior probabilities
%   out.Labels          : N×1 hard labels (argmax P)
%   out.Margin          : N×1 p_(1) - p_(2) (uncertainty gap)
%   out.Entropy         : N×1 normalized entropy in [0,1]
%   out.CriteriaTable   : table of all tried models with BIC/ICL
%   out.PlotHandles     : struct with figure handles (if plotted)
%   out.UsedVars        : variable names (columns) used from tDMC
%
% NOTES:
%   - All numeric columns are used by default. If your table has additional
%     non-DM columns, they will be ignored automatically.
%   - No scaling/normalization is performed by design (as requested).
%   - Model selection: first by max BIC; ties broken by higher ICL.
%
% Author: (you)
% -------------------------------------------------------------------------

%% Parse inputs
p = inputParser;
p.addParameter('KRange', 1:10);
p.addParameter('CovarianceTypes', {'full','diag'});
p.addParameter('SharedCovariance', [false true]);
p.addParameter('Replicates', 10);
p.addParameter('Regularization', 1e-6);
p.addParameter('MaxIter', 1000);
p.addParameter('PlotDims', [1 2 3]);
p.addParameter('RandomState', 42);
p.parse(varargin{:});
cfg = p.Results;

if ~isempty(cfg.RandomState), rng(cfg.RandomState); end

%% Select numeric variables from tDMC
isNum = varfun(@isnumeric, tDMC, 'OutputFormat','uniform');
if ~any(isNum)
    error('tDMC must contain at least one numeric column with DM coordinates.');
end
X = table2array(tDMC(:, isNum));              % N×D
varNames = tDMC.Properties.VariableNames(isNum);
[N, D] = size(X);

% Remove rows with NaNs (GMM cannot handle NaNs)
valid = all(~isnan(X),2);
if ~all(valid)
    warning('Removed %d rows with NaNs.', sum(~valid));
    X = X(valid,:);
end

%% Grid search over K × covariance settings with multiple replicates
opts = statset('MaxIter', cfg.MaxIter, 'TolFun', 1e-6);
row = 0;
rec = []; rec(1).dummy = [];  % prealloc struct array pattern
models = cell(0);

for K = cfg.KRange
    for ct = cfg.CovarianceTypes
        for sc = cfg.SharedCovariance
            row = row + 1;
            rec(row).K               = K;           %#ok<*AGROW>
            rec(row).CovarianceType  = ct{1};
            rec(row).SharedCovariance= logical(sc);
            rec(row).BIC             = -inf;
            rec(row).ICL             = -inf;
            rec(row).nlogL           = NaN;
            rec(row).Status          = "fail";
            models{row} = [];
            try
                gm = fitgmdist(X, K, ...
                    'CovarianceType', ct{1}, ...
                    'SharedCovariance', logical(sc), ...
                    'RegularizationValue', cfg.Regularization, ...
                    'Start','plus', ...
                    'Replicates', cfg.Replicates, ...
                    'Options', opts);

                % Compute BIC and ICL (ICL = BIC + sum_ik tau_ik log tau_ik)
                [P, ~] = posterior(gm, X);
                BIC = gm.BIC;
                ICL = BIC + sum(P(:) .* log(max(P(:), eps)));

                rec(row).BIC   = BIC;
                rec(row).ICL   = ICL;
                rec(row).nlogL = gm.NegativeLogLikelihood;
                rec(row).Status= "ok";
                models{row}    = gm;
            catch ME
                warning('fitgmdist failed for K=%d, type=%s, shared=%d: %s', ...
                        K, ct{1}, logical(sc), ME.message);
            end
        end
    end
end

% Build criteria table
CriteriaTable = struct2table(rec);
% Pick best: max BIC, break ties by ICL
[~, idxBIC] = max(CriteriaTable.BIC);
candidates  = find(CriteriaTable.BIC == CriteriaTable.BIC(idxBIC));
[~, j]      = max(CriteriaTable.ICL(candidates));
bestIdx     = candidates(j);
gmBest      = models{bestIdx};
if isempty(gmBest)
    error('All fits failed. Consider reducing D, increasing Regularization, or trying diag covariance.');
end

%% Posteriors, labels, uncertainty
[P, ~] = posterior(gmBest, X);                        % N×K
[~, yHard] = max(P, [], 2);

% Uncertainty metrics
Psort   = sort(P, 2, 'descend');
Margin  = Psort(:,1) - Psort(:,2);                    % gap p1 - p2 in [0,1]
Kbest   = gmBest.NumComponents;
Entropy = -sum(P .* log(max(P, eps)), 2) / log(Kbest); % normalized 0..1

% %% Basic plots (if at least 2 dims available)
% PlotHandles = struct();
% if D >= 3 && all(cfg.PlotDims <= D)
%     dims = cfg.PlotDims(:)';
%     cBase = colorPalette(max(Kbest, 6));              % base palette
%     % Compute color per point: base hue from argmax, saturation/value from margin
%     RGB = colorByMargin(P, cBase(1:Kbest,:));
% 
%     f1 = figure('Name','DMC: clusters and uncertainty','Color','w');
%     scatter3(X(:,dims(1)), X(:,dims(2)), X(:,dims(3)), 24, RGB, 'filled'); grid on; box on
%     xlabel(varNames{dims(1)}, 'Interpreter','none');
%     ylabel(varNames{dims(2)}, 'Interpreter','none');
%     zlabel(varNames{dims(3)}, 'Interpreter','none');
%     axis tight
%     title(sprintf('GMM: K=%d, %s, sharedCov=%d', ...
%           gmBest.NumComponents, gmBest.CovarianceType, gmBest.SharedCovariance));
% 
%     hold on
%     try
%         % Draw 50% and 95% Gaussian ellipses for visual shape (on selected dims)
%         for k = 1:Kbest
%             mu = gmBest.mu(k, dims);
%             if strcmp(gmBest.CovarianceType,'diag')
%                 S  = diag(gmBest.Sigma(dims, k).^2); % diag stores std? gmdistribution stores variances, not std.
%                 % Correct: for 'diag', Sigma is D×K with variances per dim
%                 S  = diag(gmBest.Sigma(dims, k));    % variances along dims
%             else
%                 % 'full': Sigma is D×D×K
%                 S  = gmBest.Sigma(dims, dims, k);
%             end
%             plot_gauss2d(mu, S, 0.50, cBase(k,:), 2.0);
%             plot_gauss2d(mu, S, 0.95, cBase(k,:)*0.7, 1.0);
%         end
%     catch
%         % If Sigma unpacking fails, silently skip ellipses
%     end
%     hold off
%     PlotHandles.Scatter = f1;
% 
%     % BIC/ICL vs K
%     f2 = figure('Name','Model selection: BIC/ICL','Color','w');
%     tiledlayout(2,1,'TileSpacing','compact');
%     nexttile;
%     hold on
%     legends = {};
%     styles  = unique(CriteriaTable(:,{'CovarianceType','SharedCovariance'}));
%     for i = 1:height(styles)
%         mask = strcmp(CriteriaTable.CovarianceType, styles.CovarianceType{i}) & ...
%                CriteriaTable.SharedCovariance == styles.SharedCovariance(i) & ...
%                CriteriaTable.Status == "ok";
%         if ~any(mask), continue; end
%         [Ks, ord] = sort(CriteriaTable.K(mask));
%         plot(Ks, CriteriaTable.BIC(mask)); %#ok<*NOPRT>
%     end
%     % (To avoid long code here, we’ll plot simpler:)
%     clf(f2);
%     f2 = figure('Name','Model selection: BIC/ICL','Color','w');
%     subplot(2,1,1);
%     plot(CriteriaTable.K, CriteriaTable.BIC, 'o-'); grid on; xlabel('K'); ylabel('BIC');
%     title('BIC across models');
%     subplot(2,1,2);
%     plot(CriteriaTable.K, CriteriaTable.ICL, 'o-'); grid on; xlabel('K'); ylabel('ICL');
%     title('ICL across models');
%     PlotHandles.Criteria = f2;
% end

%% Pack outputs
out = struct();
out.gm            = gmBest;
out.P             = P;
out.labels        = yHard;
out.K             = Kbest;
out.Margin        = Margin;
out.Entropy       = Entropy;
out.CriteriaTable = CriteriaTable;
out.UsedVars      = varNames;

end

%% ---------- Helpers (local functions) -----------------------------------




function plot_gauss2d(mu,S,p,color,lw)
%PLOT_GAUSS2D  Draw constant-probability ellipse for 2D Gaussian.
%   mu : 1×2 mean
%   S  : 2×2 covariance matrix (variances on diagonal)
%   p  : enclosed probability mass (e.g., 0.50 or 0.95)
%   color : 1×3 RGB
%   lw    : line width
if nargin<5, lw=1.5; end
r = sqrt(chi2inv(p,2));
[V,D] = eig(S);
t = linspace(0,2*pi,200);
u = [cos(t); sin(t)];
E = (V*sqrt(D)) * (r*u);
E = E' + mu;
plot(E(:,1), E(:,2), '-', 'Color', color, 'LineWidth', lw);
end
