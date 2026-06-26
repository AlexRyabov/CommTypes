function C = colorPalette(K)
%COLORPALETTE  Colorblind-friendly base palette (Tableau-like).
base = [0.121 0.466 0.705;   % blue
        1.000 0.498 0.054;   % orange
        0.172 0.627 0.172;   % green
        0.839 0.153 0.157;   % red
        0.580 0.404 0.741;   % purple
        0.549 0.337 0.294;   % brown
        0.890 0.467 0.761;   % pink
        0.498 0.498 0.498;   % gray
        0.737 0.741 0.133;   % olive
        0.090 0.745 0.811];  % cyan
if K <= size(base,1)
    C = base(1:K,:);
else
    % If more colors are needed, cycle and jitter slightly
    t = linspace(0,1,K)';
    C = hsv2rgb([t, 0.6*ones(K,1), 0.9*ones(K,1)]);
end
end