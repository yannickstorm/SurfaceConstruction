function FullKder = ComputeCovMatFullTP(R,locations,noiseVals,noiseGrad)

[D,N] = size(locations);
FullKder = CovMatFullTP(R,locations);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));