function covMatStarGrad = covMatStarGrad(sigma,gamma,X1,X2)
%% Description

%% Inputs
% sigma - GP stdev
% gamma = 1/L^2 (L - GP lengthscale)
% X1 - input locations 1
% X2 - input locations 2

%% Outputs
% KderX1X2 - Covariance Matrix


[D,N1] = size(X1);
N2 = size(X2,2);
cov = @(x1,x2)...
    (sigma^2 * exp(-1/2 * gamma *(x1 - x2)'*(x1 - x2)));

covMatStarGrad = zeros(N1 * D, N2 * (D + 1));

for n1 = 1 : N1
    for d1 = 1 : D
        for n2 = 1 : N2
            for d2 = 0 : D
                covMatStarGrad((n1 - 1) * D + d1,...
                    (n2 - 1) * (D + 1) + d2 + 1) = ...
                    covxixj(gamma,X1(:,n1),X2(:,n2),d1,d2,cov);
            end
        end
    end
end


function covxixj = covxixj(gamma,x1,x2,i,j,cov)
if (i == 0) && (j == 0)
    covxixj = cov(x1,x2);
elseif (i > 0) && (j == 0)
    covxixj = - gamma * (x1(i)-x2(i)) * cov(x1,x2);
elseif (i == 0) && (j > 0)
    covxixj = gamma * (x1(j)-x2(j)) * cov(x1,x2);
else
    covxixj = (gamma * (i==j) - ...
    gamma^2 * (x1(i)-x2(i)) * (x1(j)-x2(j))) * ...
    cov(x1,x2);
end