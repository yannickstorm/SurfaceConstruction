function Fplus = ComputeFplus(locations,surfNormals,meanValue,meanGrad, WN)

[D,~] = size(locations);
D = D * WN;

MeanGrad = ComputeMeanGrad(locations,meanValue,meanGrad, WN);

if WN
    MeanGradData(2:D+1,:) = surfNormals;
    Fplus = MeanGradData(:) - MeanGrad;
else
    Fplus = - MeanGrad;
end


function Fplus = ComputeMeanGrad(locations,meanValue,meanGrad, WN)

[D,N] = size(locations);
D = D * WN;

Fplus = zeros(D + 1,N);
for n = 1:N
    Fplus(1,n) = meanValue(locations(:,n));
    if WN
        Fplus(2:D+1,n) = meanGrad(locations(:,n));
    end
end
Fplus = Fplus(:);