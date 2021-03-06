close all
addpath('covMat')

noiseVals = 0.00000;
noiseGrad = 0.037;
sigma = 0.1243;%0.2844; %0.1
L0 = .2823; 
gamma = 2.5977;%0.6525;%1/L0^2;

% % load '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_total'
% load  '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_11'
load 'bmw_total'
% load 'bmw_11'
% load  'bmw_11'
Prior.kernel = 'SE';
locations = PartMeans;
surfNormals = SurfNormals;

% % For bmw_11
% cx = -6.1655;
% cy = -0.0472;
% cz = -3.6693;
% 
% a = 1.0763;
% b = 2.3691;
% c = 1.3254;
% r = [2 * pi , 0, 0];
% prior_type = 'E'
% 
% Prior = struct('pos',[cx cy cz],'type',prior_type, 'param', [a b c], 'rot', r);

[meanValue, meanGrad] = computePriorFunctions(Prior)

dist = 0.3;
initPoints = locations;%r * [-4;0;-2];

JMax = 1000;

[faces, vertices] = computeSurface(PartMeans, SurfNormals, ...
    Prior, ...
    meanValue, meanGrad, initPoints, dist, JMax, false);



figure
hold on
set(gca,'dataaspectratio',[1 1 1])

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');

cropCar = false;
if cropCar
    thresh = -2;
    zlim([thresh;max(vertices(:,3)) + 1]);
end

patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','gouraud');
camlight('headlight')
set(gca,'view',[46.8000   18.8000]);
% light('Position',[-1 -1 0])
