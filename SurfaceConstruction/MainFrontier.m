function [faces, vertices] = computeSurface(Locations, Gradients, covFun, initPoint, meanValue, meanGradient)

fullCovMat = computeFullCovMat(Locations, Gradients, covFun);
covMatStarValue = computeCovMatStarValue(Locations, Gradients, covFun, x);
covMatStarGrad = computeCovMatStarGrad(Locations, Gradients, covFun, x);
fGP = @(x)()


f = @(x)(meanValue(x) + fGP(x));
grad = @(x)(meanGrad(x) + gradGP(x));

% f = @(x)([1 0 0] * x);
% fd = @(x)([1;0;0]);


dist = 0.8;
figure
axis equal
hold on
x(:,1) = [r; 0; 0];
x(:,2) = x(:,1) + [0; dist; 0];
x(:,2) = Newton(x(:,2), f, grad);
x(:,3) = thirdPoint(x(:,1), x(:,2), ...
    grad(x(:,1)), grad(x(:,1)), ...
    dist, 1);
x(:,3) = Newton(x(:,3), f, grad);

plot3(x(1,[1:3 1]) , x(2,[1:3 1]), x(3,[1:3 1]), 'bo-')
plot3(0,0,0,'go')

frontiers{1}.inds = [1 2 3];
frontiers{1}.numPts = 3;
numFrontiers = 1;
nMax = 200;
faces = [1 2 3];
j = 1;
numPts = 3;
set(gca,'view',[-56.4000   24.4000]);
removeFrontiers = [];
newFrontiers = {};
numNewFrontiers = 0;
while numFrontiers > 0 && j < nMax
    
    for k = 1:numFrontiers
        index1 = frontiers{k}.inds(end);
        index2 = frontiers{k}.inds(1);
        
        xCand = thirdPoint(...
            x(:,index1),x(:,index2), ...
            grad(x(:,index1)), grad(x(:,index2)), ...
            dist, -1);
        nearIndex = check(xCand, x(:,frontiers{k}.inds), 0.8 * dist);
        if (nearIndex == 0)
            newIndex = numPts + 1;
            faces = [faces; [index1, newIndex, index2]];
            x(:,newIndex) = Newton(xCand, f, grad);
            frontiers{k}.inds = [frontiers{k}.inds newIndex];
            plot3(x(1,[index1 newIndex index2]), ...
                x(2,[index1 newIndex index2]), ...
                x(3,[index1 newIndex index2]), 'bo-');
            numPts = newIndex;
            frontiers{k}.numPts = frontiers{k}.numPts + 1;
        elseif nearIndex == 2
            faces = [faces; [index1, frontiers{k}.inds(2), index2]];
            plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                x(2,[index1 frontiers{k}.inds(2)]), ...
                x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
            
            frontiers{k}.inds(1) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            
        else
            numNewFrontiers = numNewFrontiers + 1;
            faces = [faces; [index1, frontiers{k}.inds(nearIndex), index2]];
            plot3(x(1,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                x(2,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                x(3,[index1 frontiers{k}.inds(nearIndex) index2]), 'bo-');
            newFrontiers{numNewFrontiers}.inds = frontiers{k}.inds(nearIndex:end);
            newFrontiers{numNewFrontiers}.numPts = length(newFrontiers{numNewFrontiers}.inds);
            
            frontiers{k}.inds = frontiers{k}.inds(1:nearIndex);
            frontiers{k}.numPts = length(frontiers{k}.inds);
            
        end
        
    end
    frontiers = [frontiers newFrontiers];
    numFrontiers = numFrontiers + numNewFrontiers;
    newFrontiers = {};
    numNewFrontiers = 0;
    
    %%%
    for k = 1:numFrontiers
        frontiers{k}.plot = plotFrontier(gca, frontiers{k}, x);
    end
    drawnow
    for k = 1:numFrontiers
        delete(frontiers{k}.plot);
    end
    %%%
    
    j = j + 1;
    
    
    for k = 1:numFrontiers
        if frontiers{k}.numPts < 3
            removeFrontiers = [removeFrontiers k];
        elseif frontiers{k}.numPts == 3
            removeFrontiers = [removeFrontiers k];
            faces = [faces; frontiers{k}.inds];
        end
        
    end
    frontiers(removeFrontiers) = [];
    numFrontiers = numFrontiers - length(removeFrontiers);
    removeFrontiers = [];
end
vertices = x';
patch('faces',faces,'vertices',vertices,'facecolor','green', 'edgecolor', 'none');
% 
% e(2:3,:) = [Index1, IndexN; Index2, IndexN];
% 
% plot3(x(1,IndexN), x(2,IndexN), x(3,IndexN), 'r.')
