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

covMatFull = zeros(N1 * (D + 1),N1 * (D + 1));

for n1 = 1 : N1
    for n2 = 1 : N1
        dist = X(:,n1) - X(:,n2);
        covx1x2 = sigma^2 * exp(-1/2 * gamma *...
            (dist'*dist));
        for d1 = 0 : D
            for d2 = 0 : D
                if ((n1 - 1) * (D + 1) + d1 + 1 <= ...
                        (n2 - 1) * (D + 1) + d2 + 1)
                    if (d1 == 0) && (d2 == 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            covx1x2;
                    elseif (d1 > 0) && (d2 == 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            - gamma * dist(d1) * covx1x2;
                    elseif (d1 == 0) && (d2 > 0)
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            gamma * dist(d2) * covx1x2;
                    else
                        covMatFull((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            (gamma * (d1==d2) - ...
                            gamma^2 * dist(d1) * dist(d2)) * ...
                            covx1x2;
                    end
                end
            end
        end
    end
end
covMatFull = covMatFull + covMatFull' - diag(diag(covMatFull));

