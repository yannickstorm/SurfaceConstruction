function p = ComputeGPProb(Object,Part,derIndexn)

Query = Object.Fplus(derIndexn);
% 
% if any(abs(Query) > 5 * Part.Threshold)
%    p = 0;
%     return
% end 


P = Object.P(:,derIndexn);
Q = Object.Q;

FstarMean = P' * Q;


D = numel(Query);

Tempa = P' * P;
Temp1 = Part.Kder_starstar - Tempa;
Rchol = chol(Temp1);
Det = det(Rchol)^2;
Temp = Rchol \ (Query - FstarMean);
a = Temp' * Temp;

% FstarCovMat = Part.Kder_starstar - P' * P;
% Det = det(FstarCovMat);
% Temp = (Query - FstarMean);
% a = Temp' * (FstarCovMat \ Temp);

p = exp(-1/2 * (a));
p = p * (2*pi)^(-D/2) * Det^(-1/2);


if Object.Ind == 1
    
end