%Clean variables, close windows, initialize all project constants and parameters
%I choose radial clusters because they to no give so much elongated
%communities

proj_Initialize;
tDMC = readtable(fullfile(HGP.Data.Path2ResultsData, 'Samples_DMC.xlsx'));
dmdim2use = 1:6;

gmm_out = gmm_cluster(tDMC(:, dmdim2use), 'KRange', 9:10, 'Regularization', 1e-6); %find optimal gaussian clustering

%% Save Best gaussian model 
model = gmm_out.gm;                          % your best model
meta.UsedVars = struct('Description', 'best_gmm', 'dmdim2use', dmdim2use);         % column names used to fit
meta.Date     = datestr(now, 30);     % ISO-ish timestamp
%save('cluster_best_model.mat', 'gm', 'meta');
save(fullfile(HGP.Data.Path2ResultsData, 'clust_model.mat'), 'model', 'meta');





%% Posteriors, labels, uncertainty
[P, ~] = posterior(model, tDMC{:, dmdim2use});                        % N×K
[~, yHard] = max(P, [], 2);

% Uncertainty metrics
Psort   = sort(P, 2, 'descend');
Margin  = Psort(:,1) - Psort(:,2);                    % gap p1 - p2 in [0,1]
Kbest   = model.NumComponents;
Entropy = -sum(P .* log(max(P, eps)), 2) / log(Kbest); % normalized 0..1
% Compute color per point: base hue from argmax, saturation/value from margin
cBase = f_Clrs_fresh(Kbest);  
RGB = colorByMargin(P, cBase(1:Kbest,:));

f1 = figure('Name','DMC: clusters and uncertainty','Color','w');
clf
scatter3(tDMC{:,1}, tDMC{:,2}, tDMC{:,3},...
    40, RGB, 'filled', 'MarkerEdgeColor', 'k');


%%
try
    % Draw 50% and 95% Gaussian ellipses for visual shape (on selected dims)
    dims = [1, 2];
    % dims = [1 2];  % the two coordinates you want to plot
    ctype  = lower(model.CovarianceType);           % 'diagonal' or 'full'
    shared = model.SharedCovariance;                % true/false

    for k = 1:model.NumComponents
        mu = model.mu(k, dims);

        if startsWith(ctype,'diag') || strcmp(ctype,'diagonal')
            if shared
                % Sigma is a 1×D or D×1 vector of variances (shared across components)
                v = model.Sigma(:);                 % force column D×1
                S = diag(v(dims));                   % 2×2 covariance for dims
            else
                % Sigma is D×K (each column k = variances of component k)
                v = model.Sigma(:, k);              % D×1
                S = diag(v(dims));                   % 2×2
            end
        else % 'full'
            if shared
                % Sigma is D×D (shared across components)
                S = model.Sigma(dims, dims);        % 2×2
            else
                % Sigma is D×D×K
                S = model.Sigma(dims, dims, k);     % 2×2
            end
        end

        % draw 50% and 95% ellipses
        hold on
        plot_gauss2d(mu, S, 0.50, cBase(k,:), 2.0);
        plot_gauss2d(mu, S, 0.95, cBase(k,:)*0.7, 1.0);
        hold off
    end

catch
    % If Sigma unpacking fails, silently skip ellipses
end
hold off



figure(2)
clf
nexttile

clrs = f_Clrs_fresh(gmm_out.gm.NumComponents)
%colorPalette(out_rad_cl.K);

hold on
for k = 1:gmm_out.K
    idx = (gmm_out.labels == k);
    scatter3(tDMC{idx,1}, tDMC{idx,2}, tDMC{idx,3}, ...
        40, repmat(clrs(k,:), sum(idx),1), 'filled', 'MarkerEdgeColor', 'k');
end
legend('Location','best')
hold off
title(['Gaussian mixture with ' NS(gmm_out.gm.NumComponents) ' clusters'])
%%
out_rad_cl = radius_cluster(tDMC{:, 1:6}, 'Radius', 0.5, 'Manner', 'approx', 'ReturnCenters', true);
clrs = f_Clrs_fresh(out_rad_cl.K)

nexttile(2)
cla
hold on
for k = 1:out_rad_cl.K
    idx = (out_rad_cl.labels == k);
    scatter3(tDMC{idx,1}, tDMC{idx,2}, tDMC{idx,3}, ...
        40, repmat(clrs(k,:), sum(idx),1), 'filled', 'MarkerEdgeColor', 'k');
end
legend('Location','best')
hold off
% 4) check how compact are the clusters
fprintf('K=%d, median rad90=%.3g, sizes in [%d, %d]\n', ...
    out_rad_cl.K, median(out_rad_cl.rad90), min(out_rad_cl.sizes), max(out_rad_cl.sizes));

title(sprintf('Radius clustering (R=%.3g, K=%d)', out_rad_cl.radius, out_rad_cl.K));

%%

