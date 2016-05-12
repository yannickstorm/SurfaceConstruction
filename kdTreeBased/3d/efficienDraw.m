% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

load('bmw_11.mat');

X = PartMeans;
m = length(X); 
centerData = sum(PartMeans')'/m;
for i =1:m
   X(:,i) = X(:,i) - centerData; 
end

data = zeros(1,length(SurfNormals));
data = [data ; SurfNormals]; 
f = reshape(data, [1,4*length(X)])';

plot3(X(1,:), X(2,:),X(3,:), 'r.', 'MarkerSize',40);
hold on;
quiver3(X(1,:), X(2,:), X(3,:), data(2,:), data(3,:),data(4,:));


sigma = 0.12;
gamma = 2.6;
noiseVal = 0.001;
noiseGrad = 0.001;
display('Computing means');
R = 1;
cen = [0,0,0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);
meandx = @(x) 1/R*((x(1)-cen(1)));
meandy = @(x) 1/R*((x(2)-cen(2)));
meandz = @(x) 1/R*((x(3)-cen(3)));

display('Computing mean and covariance of data');
K = ComputeFullKder(sigma, gamma, X, noiseVal, noiseGrad);
for i = 1:m
    mu((i-1)*4 +1) = mean(X(:,i));
    mu((i-1)*4 +2) = meandx(X(:,i));
    mu((i-1)*4 +3) = meandy(X(:,i));
    mu((i-1)*4 +4) = meandz(X(:,i));
end

%% Efficiend draw
iterations = 3;

xLimits = [-2, 2];
yLimits = [-2, 2];
zLimits = [-2, 2];
centroid = [sum(xLimits)/2, sum(yLimits)/2,sum(zLimits)/2];

% 2 valid, 1 new, 0 invalid.
root = {2,  {}, xLimits, yLimits, centroid, -1, zLimits};
Qmat = inv(K)*(f - mu');

evalFun = @(x) [mean(x), meandx(x), meandy(x), meandz(x)]' + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;
tic
root = expandCell(root,evalFun, 2);
root = expandCell(root,evalFun, 2);
root = expandCell(root,evalFun, 2);
for iter=1:iterations
    % Expand tree
    root = expandCell(root, evalFun, 1); 
    % Validate branches
    root = validatePoints(root, root, evalFun);
    
    figure();points = [];
    cols = [];
    plot3(X(1,:), X(2,:),X(3,:), 'r.', 'MarkerSize',40);
    hold on;
    quiver3(X(1,:), X(2,:), X(3,:), data(2,:), data(3,:),data(4,:));
    [points, cols] = getPointsTree(points, cols, root,false);
    plot3(points(:,1), points(:,2), cols, 'o')
    pause();
    close all;
    
    
end
root = validatePoints(root, root, evalFun);
toc
figure();
hold on;
points = [];
cols = [];
surface = getPointsTree(points, cols, root, true);
plot3(surface(:,1), surface(:,2),surface(:,3), 'go');
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
quiver3(X(1,:), X(2,:), X(3,:), data(2,:), data(3,:),data(4,:));
grid;
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
axis equal


% Ground Truth
[Xg,Yg, Zg] = meshgrid(-2:0.5:2,-1.4:0.5:2,-2:0.5:2);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1),reshape(Zg,d1*d2,1)]';
n = length(Xs);

for i = 1:n
    mus((i-1)*4 +1) = mean(Xs(:,i));
    mus((i-1)*4 +2) = meandx(Xs(:,i));
    mus((i-1)*4 +3) = meandy(Xs(:,i));
    mus((i-1)*4 +4) = meandz(Xs(:,i));
end
mus = mus';
Ks = ComputeKderX1X2(sigma, gamma, Xs, X)';
Kss = ComputeFullKder(sigma, gamma, Xs, noiseVal, noiseGrad);

fs = mus + Ks'*inv(K)*(f - mu');
sig = Kss' - Ks'*inv(K)*Ks;

Fs = reshape(fs(1:4:end),7,9,9);

figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
figure();
hold on;
pMean = patch(  isosurface(Xg, Yg, Zg, Fs, 0), ...
                'FaceColor','green',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
quiver3(X(1,:), X(2,:), X(3,:), data(2,:), data(3,:),data(4,:));