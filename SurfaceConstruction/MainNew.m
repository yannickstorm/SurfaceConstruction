close all
addpath('covMat')

noiseVals = 0.00000;
noiseGrad = 0.037;
sigma = 0.1243;%0.2844; %0.1
L0 = .2823; 
gamma = 2.5977;%0.6525;%1/L0^2;

% % load '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_total'
% load  '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_11'
% load 'bmw_total'
load  'bmw_11'
% For bmw_11
cx = -6.1655;
cy = -0.0472;
cz = -3.6693;

a = 1.0763;
b = 2.3691;
c = 1.3254;
locations = PartMeans;
surfNormals = SurfNormals;


%  
% locations = r * [[1;0;0] [0;1;0] [0;0;1]];
% surfNormals = [[0.7;0.7;0] [0;1;0] [0;1;0]];
% %Simple case
% 
% cx = 0;
% cy = 0;
% cz = 0;
% 
% a = 1;
% b = 1;
% c = inf;

A = diag([1/a^2, 1/b^2, 1/c^2]);
loc = [cx cy cz]';
meanValue = @(x)(a/2 * ((x-loc)' * A * (x-loc) - 1));
meanGrad = @(x)(a * A * (x-loc));

dist = 0.4;
initPoint = locations(:,1);%r * [-4;0;-2];

figure
hold on
axis equal

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');



[faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist, false);

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
set(gca,'view',[46.8000   18.8000]);