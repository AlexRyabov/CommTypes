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