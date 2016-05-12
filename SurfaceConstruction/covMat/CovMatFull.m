function covMatFull = CovMatFull(sigma,gamma,X1)
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


[D,N1] = size(X1);
cov = @(x1,x2)...
    (sigma^2 * exp(-1/2 * gamma *(x1 - x2)'*(x1 - x2)));

covMatFull = zeros(N1 * (D + 1),N1 * (D + 1));

for n1 = 1 : N1
    for n2 = 1 : N1
        covx1x2 = cov(X1(:,n1),X1(:,n2));
        for d1 = 0 : D
            for d2 = 0 : D
                if ((n1 - 1) * (D + 1) + d1 + 1 <= ...
                    (n2 - 1) * (D + 1) + d2 + 1)
                    covMatFull((n1 - 1) * (D + 1) + d1 + 1,...
                        (n2 - 1) * (D + 1) + d2 + 1) = ...
                        covxixj(gamma,X1(:,n1),X1(:,n2),covx1x2,d1,d2);
                end
            end
        end
    end
end
if isequal(X1,X1)
    covMatFull = covMatFull + covMatFull' - diag(diag(covMatFull));
end

function covxixj = covxixj(gamma,x1,x2,covx1x2,i,j)
if (i == 0) && (j == 0)
    covxixj = covx1x2;
elseif (i > 0) && (j == 0)
    covxixj = - gamma * (x1(i)-x2(i)) * covx1x2;
elseif (i == 0) && (j > 0)
    covxixj = gamma * (x1(j)-x2(j)) * covx1x2;
else
    covxixj = (gamma * (i==j) - ...
    gamma^2 * (x1(i)-x2(i)) * (x1(j)-x2(j))) * ...
    covx1x2;
end