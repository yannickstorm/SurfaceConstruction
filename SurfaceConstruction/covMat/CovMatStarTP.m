function covMatStar = CovMatStarTP(R,X1,X2)
%% Description
% Computes the covariance matrix between function values and 1st order
% derivatives between input locations X1 and X2
% given GP parameters sigma, gamma
% X1 = X2 returns the covariance matrix between function values and 1st order
% derivatives at X1

%% Inputs
% R
% X1 - input locations 1
% X2 - input locations 2

%% Outputs
% KderX1X2 - Covariance Matrix


[D,N1] = size(X1);
N2 = size(X2,2);
covMatStar = zeros(N1 * (D + 1),N2 * (D + 1));

for n1 = 1 : N1
    for n2 = 1 : N2
        dist = X1(:,n1) - X2(:,n2);
        distNorm = norm(dist);
        for d1 = 0 : D
            for d2 = 0 : D
                if (d1 == 0) && (d2 == 0)
                    covMatStar((n1 - 1) * (D + 1) + d1 + 1,...
                        (n2 - 1) * (D + 1) + d2 + 1) = ...
                        2 * distNorm^3 ...
                        - 3 * R * distNorm^2 + R^3;
                elseif (d1 > 0) && (d2 == 0)
                    covMatStar((n1 - 1) * (D + 1) + d1 + 1,...
                        (n2 - 1) * (D + 1) + d2 + 1) = ...
                        6 * dist(d1) * (distNorm - R);
                elseif (d1 == 0) && (d2 > 0)
                    covMatStar((n1 - 1) * (D + 1) + d1 + 1,...
                        (n2 - 1) * (D + 1) + d2 + 1) = ...
                        - 6 * dist(d2) * (distNorm - R);
                else
                    if X1(:,n1) == X2(:,n2)
                        covMatStar((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            - 6 * ((d1 == d2) * (distNorm - R));
                    else
                        covMatStar((n1 - 1) * (D + 1) + d1 + 1,...
                            (n2 - 1) * (D + 1) + d2 + 1) = ...
                            - 6 * (dist(d1) * dist(d2)/distNorm ...
                            + (d1 == d2) * (distNorm - R));
                    end
                end
            end
        end
    end
end
