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
load 'eggplant2'
% load  'bmw_11'
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

dist = 0.4;
initPoints = locations(:,42);%r * [-4;0;-2];


[faces, vertices, frontiers] = computeSurface(locations, surfNormals, ...
    Prior, ...
    meanValue, meanGrad, initPoints, dist, false, 485);



figure
hold on
set(gca,'dataaspectratio',[1 1 1])
% 
% plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
% quiver3(locations(1,:),locations(2,:),locations(3,:),...
%     surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');
% quiver3(locations(1,:),locations(2,:),locations(3,:),...
%     surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');
% 


patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', [0.8 0.8 0.8], ...
    'facelighting','gouraud');
camlight('headlight')
set(gca,'view',[-2.0000   38.0000], 'visible', 'off');
for k = 1:length(frontiers)
    plotFrontier(gca, frontiers{k}, vertices');
end

% light('Position',[-1 -1 0])
