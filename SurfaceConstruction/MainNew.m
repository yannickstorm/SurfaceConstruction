close all
addpath('covMat')

sigma =   0.0058;
gamma =  688.7485;
noiseVals =   0.000001;
noiseGrad =    0.0001;
r = 1;

locations = r * [[2;0;0] [0;1;0]];
surfNormals = [[1;0;0] [0;1;0]];

A = 1/r^2 * eye(3);

meanValue = @(x)(r/2 * (x' * A * x - 1));
meanGrad = @(x)(r * A * x);

dist = 0.4;
initPoint = r * [1;0;0];

[faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist);