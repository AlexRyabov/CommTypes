function out = radius_cluster(X, varargin)
%RADIUS_CLUSTER_DMC  Call Python 'radius-clustering' on DM coordinates.
% INPUT:
%   X : N×D double matrix (Diffusion Map coordinates; no scaling inside).
%
% NAME-VALUE (all optional):
%   'Radius'        : positive scalar, required for the algorithm (e.g., 0.15)
%   'Manner'        : 'approx' (default) or 'exact'
%   'Metric'        : distance name (currently Euclidean only in the Python lib’s standard API).
%   'RandomState'   : integer seed (optional; passed when supported)
%   'ReturnCenters' : true/false, compute cluster medoids in MATLAB (default: true)
%
% OUTPUT struct:
%   out.labels      : N×1 double cluster labels (1..K)
%   out.K           : number of clusters
%   out.radius      : radius used
%   out.manner      : 'approx'/'exact'
%   out.centers     : K×D medoids (if ReturnCenters=true)
%   out.sizes       : K×1 cluster sizes
%   out.rad90       : K×1 90th percentile distance-to-center (post-check)
%
% Notes:
%   - Requires Python 'radius-clustering' and 'numpy' available in pyenv.
%   - The Python package labels are 0-based; we convert to 1-based.
%   - If you need a custom metric, pre-transform X or use DM space as-is.

% Parse args
p = inputParser;
p.addParameter('Radius', [], @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Manner', 'approx', @(s)ischar(s)||isstring(s));
p.addParameter('Metric', 'euclidean', @(s)ischar(s)||isstring(s));
p.addParameter('RandomState', [], @(x)isempty(x)||(isscalar(x)&&isnumeric(x)));
p.addParameter('ReturnCenters', true, @(x)islogical(x)||ismember(x,[0,1]));
p.parse(varargin{:});
cfg = p.Results;

if isempty(cfg.Radius)
    error('Please provide ''Radius'', e.g., radius_cluster_dmc(X, ''Radius'', 0.15)');
end
if ~strcmpi(cfg.Metric,'euclidean')
    warning('Only Euclidean is directly supported here. Using Euclidean.');
end

% Check Python availability
try
    pynp = py.importlib.import_module('numpy');
    pyrc = py.importlib.import_module('radius_clustering');
catch ME
    error(['Python packages not found. Run setup_radius_clustering.m first. ', ...
           'Original error: %s'], ME.message);
end

% Convert X -> numpy array (float64)
Xnp = pynp.array(X);

% Instantiate and fit the model
% API: RadiusClustering(manner='approx'/'exact', radius=..., **kwargs)
kwargs = pyargs('manner', string(cfg.Manner), 'radius', cfg.Radius);
rc = pyrc.RadiusClustering(kwargs);
rc.fit(Xnp);

% Fetch labels (0-based → 1-based)
lbl_py = py.getattr(rc,'labels_');
lbl = double(pynp.asarray(lbl_py));
lbl = lbl(:) + 1;

% Pack basic outputs
out = struct();
out.labels = lbl;
out.K = numel(unique(lbl));
out.radius = cfg.Radius;
out.manner = char(cfg.Manner);

% Sizes
out.sizes = accumarray(lbl, 1, [out.K, 1]);

% Optional: compute simple medoids and compactness diagnostics in MATLAB
if cfg.ReturnCenters
    Dfun = @(A,b) sqrt(sum((A - b).^2,2)); % Euclidean
    centers = zeros(out.K, size(X,2));
    rad90 = zeros(out.K,1);
    for k=1:out.K
        idx = (lbl==k);
        Xk = X(idx,:);
        % medoid: pick point with minimal sum of distances
        % (O(n^2) per cluster; OK for moderate sizes)
        if size(Xk,1) <= 2000
            D = pdist2(Xk, Xk);
            [~,midx] = min(sum(D,2));
            centers(k,:) = Xk(midx,:);
        else
            % fallback: geometric median approx via the mean
            centers(k,:) = mean(Xk,1);
        end
        d = Dfun(Xk, centers(k,:));
        rad90(k) = prctile(d, 90);
    end
    out.centers = centers;
    out.rad90 = rad90;
end
end
