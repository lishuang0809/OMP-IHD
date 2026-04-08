function [X, Q, S_opt] = data_generation(N, L, k, J)

% This function is used to generate a dictionary Q of size NxL, 
% and J k-sparse vectors (stored in S_opt), 
% and J measurement vectors (stored in X).
% Coded by Shuang Li, Feburary 2024.
% Updated by adding ||X(:,j)||_2 = ||Q S_opt(:,j)||_2 = 1, January 2025.

Q = randn(N,L);
Q = normc(Q);      % normalize the column of Q


S_opt = zeros(L, J);  % Initialize the sparse matrix

% Generate k-sparse vectors
for j = 1:J    
    loc = randperm(L, k); % Randomly select k unique positions in the column
    S_opt(loc, j) = randn(k, 1); % Assign random values to these positions
end

X = Q * S_opt;   % Generate measurement vectors
S_opt = S_opt*diag(1./vecnorm(X));
X = Q * S_opt;   % Such S_opt makes  ||X(:,j)||_2 = 1, for all j







