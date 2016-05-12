function FullKder = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad)

[D,N] = size(locations);
FullKder = CovMatFull(sigma,gamma,locations);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));