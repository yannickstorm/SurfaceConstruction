close all
addpath('covMat')

sigma =   .1;
gamma =  0.8;
noiseVals =   0.00001;
noiseGrad =    0.00001;
r = 1;

locations = r * [[1;0;0] [0;1;0] [0;0;1] [0.7;0.7;0]];
surfNormals = [[0.7;0.7;0] [0;1;0] [0;-1;0] [1;0;0]];


A = 1/r^2 * eye(3);

meanValue = @(x)(r/2 * (x' * A * x - 1));
meanGrad = @(x)(r * A * x);

dist = 0.25;
initPoint = r * [1;0;0];

[faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist);

figure
hold on
axis equal

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');


patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none');
camlight