function Fplus = ComputeFplus(locations,surfNormals,meanValue,meanGrad)

[D,~] = size(locations);
MeanGrad = ComputeMeanGrad(locations,meanValue,meanGrad);

MeanGradData(2:D+1,:) = surfNormals;
Fplus = MeanGradData(:) - MeanGrad;

function Fplus = ComputeMeanGrad(locations,meanValue,meanGrad)

[D,N] = size(locations);

Fplus = zeros(D + 1,N);
for n = 1:N
    Fplus(1,n) = meanValue(locations(:,n));
    Fplus(2:D+1,n) = meanGrad(locations(:,n));
end
Fplus = Fplus(:);