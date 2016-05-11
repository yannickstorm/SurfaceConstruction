function FullKder = ComputeFullKder(sigma,gamma,PartMeans,noiseVals,noiseGrad,With_normals)

[D,N] = size(PartMeans);
D = D * With_normals;
FullKder = ComputeKderX1X2(sigma,gamma,PartMeans,PartMeans, With_normals);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));