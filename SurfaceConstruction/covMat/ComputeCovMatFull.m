function FullKder = ComputeCovMatFull(locations, Prior)

[D,N] = size(locations);
switch Prior.kernel
    case 'TP'
    FullKder = CovMatFullTP(Prior.param(2),locations);
    
    case 'SE'
    FullKder = CovMatFullSE(Prior.Sigma, Prior.Gamma,locations);
    
    otherwise
        error('The kernel defined in Prior does not correspond to any of the possible ones.')
end
    
NoiseDiag =  [abs(Prior.noiseVals) * ones(1,N);
              abs(Prior.noiseGrad) * ones(D,N)];
FullKder = FullKder + diag(NoiseDiag(:));