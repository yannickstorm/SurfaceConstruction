function KderX1X2 = ComputeKderX1X2(sigma,gamma,X1,X2)
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
N2 = size(X2,2);
cov = @(x1,x2)...
    (sigma^2 * exp(-1/2 * gamma *(x1 - x2)'*(x1 - x2)));

KderX1X2 = zeros(N1 * (D + 1),N2 * (D + 1));

for n1 = 1 : N1
    for d1 = 0 : D
        for n2 = 1 : N2
            for d2 = 0 : D
                if ((n1 - 1) * (D + 1) + d1 + 1 <= (n2 - 1) * (D + 1) + d2 + 1) || ~isequal(X1,X2)
                    KderX1X2((n1 - 1) * (D + 1) + d1 + 1     ,   (n2 - 1) * (D + 1) + d2 + 1) = ...
                        covxixj(gamma,X1(:,n1),X2(:,n2),d1,d2,cov);
                end
            end
        end
    end
%     n1
end
if isequal(X1,X2)
    KderX1X2 = KderX1X2 + KderX1X2' - diag(diag(KderX1X2));
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