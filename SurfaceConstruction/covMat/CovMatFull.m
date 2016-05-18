function covMatFull = CovMatFull(sigma,gamma,X)
%% Description
% Computes the covariance matrix between function values and 1st order
% derivatives between input locations X1 and X2
% given GP parameters sigma, gamma
% X1 = X2 returns the covariance matrix between function values and 1st order
% derivatives at X1

%% Inputs
% sigma - GP stdev
% gamma = 1/L^2 (L - GP lengthscale)
% X1 - input locations 1
% X2 - input locations 2

%% Outputs
% KderX1X2 - Covariance Matrix


[D,N1] = size(X);

covMatFull = zeros(N1,N1);

for n1 = 1 : N1
    for n2 = 1 : N1
        dist = X(:,n1) - X(:,n2);
        covx1x2 = sigma^2 * exp(-1/2 * gamma *...
            (dist'*dist));
        covMatFull(n1, n2) = ...
                            covx1x2;
    end
end
covMatFull = covMatFull + covMatFull' - diag(diag(covMatFull));

