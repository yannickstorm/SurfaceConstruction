function DataLogLikelihood = ...
    DataLogLikelihood(PartMeans,SurfNormals,sigma,gamma,noiseVals,noiseGrad,cx,cy,cz,a,b,c)
% Outputs the current values
[sigma,gamma,noiseVals,noiseGrad,cx,cy,cz,a,b,c]

%Compute the full covariance Matrix
FullKder = ComputeFullKder(sigma,gamma,PartMeans,noiseVals,noiseGrad);

% Estaimate the centre location
% [Mean,~] = EstimateCentreGrad(...
%                 PartMeans,...
%                 SurfNormals,...
%                 R,FullKder);
            
Object.Loc = [cx,cy,cz]';%Mean;

%Compute the Fplus vector
Fplus = ComputeFplus(PartMeans,SurfNormals,Object,cx,cy,cz,a,b,c);

%argument of the exp-function
T1 = -1/2 * Fplus' * (FullKder \ Fplus);

%log of the normalisation factor
% T2 =  - 1/2 * log(det(FullKder))
T2 =  - 1/2 * sum(log(eig(FullKder)));%log(det(FullKder))

% det(FullKder)
% eig(FullKder)
%data log likelihood
DataLogLikelihood = T1 + T2