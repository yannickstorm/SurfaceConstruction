% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

%% Data
% X = [0,-0.5;
%     -0.3,-0.1;
%     0.5,-1.0;
%      -0.5,0.5;
%      0.5,0.5;
%     1,0]';
% 
% f = [   0,-cos(20/180*pi),-sin(20/180*pi),...
%         0,-cos(45/180*pi),-sin(45/180*pi),...
%         0,-cos(120/180*pi),-sin(120/180*pi),...
%         0, -cos(45/180*pi), sin(45/180*pi),...
%         0, cos(45/180*pi), sin(45/180*pi),...
%         0,1,0]';
%     
% Data
X = ([-0.301861, -3.988338;
	 -0.301861, -2.988338;
	 -1.301861, -1.988338;
	 -1.301861, -0.988338;
	 -2.301861, 0.011662;
	 -3.301861, 0.011662;
	 -4.301861, 1.011662;
	 -1.301861, -4.988338;
	 -2.301861, -4.988338;
	 -3.301861, -5.488338;
	 -3.301861, -6.388338;
	 -2.301861, -6.988338;
	 0.698139, -6.588338;
	 1.698139, -5.188338;
	 2.698139, -3.988338;
	 2.398139, -2.988338;
	 2.698139, -1.988338;
	 2.698139, -0.288338;
	 1.998139, -0.288338;
	 1.698139, 0.511662;
	 1.298139, 1.811662;
	 1.898139, 3.311662;
	 1.898139, 4.411662;
	 1.154870, 5.217616;
	 -0.071972, 5.531043;
	 -1.218218, 5.548953;
	 -2.203273, 5.190751;
	 -3.107733, 4.868369;
	 -4.191293, 4.151966;
	 -4.648001, 3.435562;
	 -4.636122, 2.137056]')/7;

f = [   0,-cos(20/180*pi),sin(20/180*pi),...
        0,-cos(0/180*pi),sin(0/180*pi),...
        0,-cos(20/180*pi),-sin(20/180*pi),...
        0, -cos(0/180*pi), sin(0/180*pi),...
        0, cos(90/180*pi), -sin(90/180*pi),...
        0,cos(90/180*pi),-sin(90/180*pi),...
        0,cos(120/180*pi),-sin(120/180*pi),...
        0, -cos(45/180*pi), sin(45/180*pi),...
        0, cos(90/180*pi), sin(90/180*pi),...
        0,-cos(45/180*pi),-sin(45/180*pi),...
        0,cos(120/180*pi),-sin(120/180*pi),...
        0, -cos(60/180*pi), -sin(60/180*pi),...
        0, cos(45/180*pi), -sin(45/180*pi),...
        0,cos(45/180*pi),-sin(45/180*pi),...
        0,-cos(120/180*pi),-sin(120/180*pi),...
        0,cos(0/180*pi),sin(0/180*pi),...
        0,-cos(120/180*pi),-sin(120/180*pi),...
        0, cos(20/180*pi), sin(20/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi),...
        0,cos(45/180*pi),sin(45/180*pi),...
        0,cos(0/180*pi),sin(0/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi),...
        0,-cos(70/180*pi),sin(70/180*pi),...
        0,-cos(70/180*pi),sin(70/180*pi),...
        0, -cos(45/180*pi), sin(45/180*pi),...
        0, -cos(45/180*pi), sin(45/180*pi),...
        0,-cos(45/180*pi),-sin(45/180*pi),...
        0, -cos(45/180*pi), -sin(45/180*pi)]';

data = reshape(f, [3,length(X)])';
m = length(X); 

plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
hold on;
quiver(X(1,:)', X(2,:)', data(:,2), data(:,3));


sigma = 0.3;
gamma = 5;
noiseVal = 0.1;
noiseGrad = 0.1;
display('Computing means');
R = 1;
cen = [0.0, 0.0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);
meandx = @(x) 1/R*((x(1)-cen(1)));
meandy = @(x) 1/R*((x(2)-cen(2)));

display('Computing mean and covariance of data');
K = ComputeFullKder(sigma, gamma, X, noiseVal, noiseGrad);
for i = 1:m
    mu((i-1)*3 +1) = mean(X(:,i));
    mu((i-1)*3 +2) = meandx(X(:,i));
    mu((i-1)*3 +3) = meandy(X(:,i));
end
mu = mu';

%% Efficiend draw
iterations = 5;

xLimits = [-1.5, 1.5];
yLimits = [-1.5, 1.5];
centroid = [sum(xLimits)/2, sum(yLimits)/2];

% 2 valid, 1 new, 0 invalid.
root = {2,  {}, xLimits, yLimits, centroid, -1};
Qmat = inv(K)*(f - mu);

evalFun = @(x) [mean(x), meandx(x), meandy(x)]' + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;
tic
root = expandCell(root,evalFun, 2);
root = expandCell(root,evalFun, 2);
% root = expandCell(root,evalFun, 2);
% root = expandCell(root,evalFun, 2);
for iter=1:iterations
    % Validate branches
    root = validatePoints(root, root, evalFun);
    
%     figure();points = [];
%     cols = [];
%     [points, cols] = getPointsTree(points, cols, root,false);
%     plot3(points(:,1), points(:,2), cols, 'o')
%     pause();
%     close all;
    
    % Expand tree
    root = expandCell(root, evalFun, 1);  
end
root = validatePoints(root, root, evalFun);
toc
figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
points = [];
cols = [];
% [points, cols] = getPointsTree(points, cols, root,false);
% plot3(points(:,1), points(:,2), cols, 'o')
quiver(X(1,:)', X(2,:)', data(:,2), data(:,3));
grid;
axis([-1.5 1.5 -1.5 1.5]);

surface = getPointsTree(points, cols, root, true);
plot(surface(:,1), surface(:,2), 'go');
%% Ground Truth
% [Xg,Yg] = meshgrid(-1.4:0.2:1.4,-1.4:0.2:1.4);
% [d1,d2] = size(Xg);
% Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
% n = length(Xs);
% 
% for i = 1:n
%     mus((i-1)*3 +1) = mean(Xs(:,i));
%     mus((i-1)*3 +2) = meandx(Xs(:,i));
%     mus((i-1)*3 +3) = meandy(Xs(:,i));
% end
% mus = mus';
% Ks = ComputeKderX1X2(sigma, gamma, Xs, X)';
% Kss = ComputeFullKder(sigma, gamma, Xs, noiseVal, noiseGrad);
% 
% fs = mus + Ks'*inv(K)*(f - mu);
% sig = Kss' - Ks'*inv(K)*Ks;
% 
% Fs = reshape(fs(1:3:end),d1,d2);
% 
% figure();
% hold on;
% plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
% quiver(X(1,:)', X(2,:)', data(:,2), data(:,3));