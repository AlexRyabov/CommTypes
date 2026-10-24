function [d, fitmodel] = f_CorrelationDimension(D, Range)
%find correlation dimension using the matrix of distances and the Range
%where to look for the dimension

% D is NxN symmetric matrix of pairwise distances
N = size(D, 1);

% Extract the upper triangular part of D, excluding the diagonal (self-distances)
D_upper = triu(D, 1); % '1' excludes the diagonal

% Convert the upper triangular matrix to a vector of unique distances
dists = nonzeros(D_upper); % This gives a column vector of distances without zeros

% Define a range of r values (adjust the number of points and range as needed)
r_values = logspace(log10(min(dists)), log10(max(dists)), 100);

% Initialize the correlation integral array
C = zeros(length(r_values), 1);

% Total number of unique pairs
num_pairs = N * (N - 1) / 2;

for k = 1:length(r_values)
    r = r_values(k);
    % Count the number of pairs with distance less than r
    C(k) = sum(dists < r) / num_pairs;
end
%% Remove zero or infinite values before taking logarithms
idx2fit = r_values>Range(1) & r_values<Range(2);
log_C = log10(C(idx2fit));
log_r = log10(r_values(idx2fit));

% Perform linear regression to estimate the correlation dimension
fitmodel = fitlm(log_r, log_C);
d = fitmodel.Coefficients.Estimate(2); % Correlation dimension estimate
se = fitmodel.Coefficients.SE(2); % Standart error 

%Plot the results
[fg, iFgNrNew] = addFindFigure('Correlation dimension estimate', -1);
clf
loglog(r_values, C, 'o');
hold on;
plot(10.^log_r, 10.^fitmodel.predict(log_r'), '-r');
xlabel('r');
ylabel('C(r)');
title(['Correlation Dimension Estimate: D = ', NS(d, 3) '\pm' NS(se, 1) ]);
grid on;
end