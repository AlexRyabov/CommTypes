function RGB = colorByMargin(P, baseC)
%COLORBYMARGIN  Assign color by argmax hue and fade by uncertainty margin.
% Inputs:
%   P     : N×K posterior probabilities
%   baseC : K×3 base RGB colors (rows in 0..1)
[~, kmax] = max(P,[],2);
Psort = sort(P,2,'descend');
margin = Psort(:,1) - Psort(:,2);          % 0..1
N = size(P,1);
RGB = zeros(N,3);
for i=1:N
    hsv = rgb2hsv(baseC(kmax(i),:));
    hsv(2) = 0.25 + 0.75*margin(i);        % saturation increases with certainty
    hsv(3) = 0.70 + 0.30*margin(i);        % value increases with certainty
    RGB(i,:) = hsv2rgb(hsv);
end
end
