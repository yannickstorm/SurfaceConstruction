close all
addpath('covMat')

sigma =   .1;
gamma =  0.5;
noiseVals =   0.00001;
noiseGrad =    0.00001;
r = 1;

locations = r * [[1;0;0] [0;1;0] [0;0;1]];
surfNormals = [[0.7;0.7;0] [0;1;0] [0;1;0]];


A = 1/r^2 * eye(3);

meanValue = @(x)(r/2 * (x' * A * x - 1));
meanGrad = @(x)(r * A * x);

dist = 0.3;
initPoint = r * [1;0;0];

[faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist);