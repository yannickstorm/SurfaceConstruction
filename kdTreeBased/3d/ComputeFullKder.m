function FullKder = ComputeFullKder(sigma,gamma,PartMeans,noiseVals,noiseGrad)

[D,N] = size(PartMeans);
FullKder = ComputeKderX1X2(sigma,gamma,PartMeans,PartMeans);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));