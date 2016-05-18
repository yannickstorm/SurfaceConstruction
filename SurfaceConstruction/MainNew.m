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

dist = .2;
initPoints = locations;%r * [-4;0;-2];


[faces, vertices] = computeSurface(locations, surfNormals, ...
    Prior, ...
    meanValue, meanGrad, initPoints, dist, true);

figure
hold on
axis equal

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');

cropCar = true;
if cropCar
    thresh = -2;
    remFaces = [];
    vertices(vertices(:,3)<thresh,3) = thresh;
    for i = 1:size(faces,1)
        if all(vertices(faces(i,:),3) == thresh)
            remFaces = [remFaces i];
        end
    end
    faces(remFaces,:) = [];
end
        

patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','flat');
camlight
set(gca,'view',[46.8000   18.8000]);
light('Position',[-1 -1 0])
