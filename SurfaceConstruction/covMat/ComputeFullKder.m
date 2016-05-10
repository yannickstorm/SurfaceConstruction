function FullKder = ComputeFullKder(sigma,gamma,locations,noiseVals,noiseGrad)

[D,N] = size(locations);
FullKder = ComputeKderX1X2(sigma,gamma,locations,locations);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));