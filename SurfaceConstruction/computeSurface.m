function [faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist, plot)

fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);

covMatData = ComputeFullKder(sigma,gamma,locations,noiseVals,noiseGrad);
RVector = covMatData\fPlusData;

fGP = @(x)(covMatStarValue(sigma, gamma, x, locations) * RVector);
gradGP = @(x)(covMatStarGrad(sigma, gamma, x, locations) * RVector);

f = @(x)(meanValue(x) + fGP(x));
grad = @(x)(meanGrad(x) + gradGP(x));

x(:,1) = initPoint;
x(:,1) = NewtonOneStep(x(:,1), f, grad);
gradX(:,1) = grad(x(:,1));

gradInit = gradX(:,1);
perpVector = cross(gradInit,rand(3,1));
normPerpVector = perpVector/norm(perpVector);
x(:,2) = x(:,1) + dist * normPerpVector;
x(:,2) = NewtonOneStep(x(:,2), f, grad);
gradX(:,2) = grad(x(:,2));

x(:,3) = thirdPoint(x(:,1), x(:,2), ...
    grad(x(:,1)), grad(x(:,1)), ...
    dist, 1);
x(:,3) = NewtonOneStep(x(:,3), f, grad);
gradX(:,3) = grad(x(:,3));

if plot
    figure
    axis equal
    hold on
    set(gca,'view',[-116.0000   -2.8000]);
    gradients = true;
    plot3(x(1,[1:3 1]) , x(2,[1:3 1]), x(3,[1:3 1]), 'bo-')
    if gradients
        quiver3(x(1,[1:3]) , x(2,[1:3]), x(3,[1:3]),...
            gradX(1,[1:3]) , gradX(2,[1:3]), gradX(3,[1:3]),0);
    end
end

frontiers{1}.inds = [1 2 3];
frontiers{1}.numPts = 3;
numFrontiers = 1;
nMax = 4000;
faces = [1 2 3];
j = 1;
numPts = 3;
removeFrontiers = [];
newFrontiers = {};
numNewFrontiers = 0;
while numFrontiers > 0 && j < nMax
    
    for k = 1:numFrontiers
        index1 = frontiers{k}.inds(end);
        index2 = frontiers{k}.inds(1);
        
        if plot
            activeEdge = plot3(x(1,[index1 index2]), ...
                x(2,[index1 index2]), ...
                x(3,[index1 index2]), 'm-','linewidth',4);
        end
        
        xCand = thirdPoint(...
            x(:,index1),x(:,index2), ...
            gradX(:,index1), gradX(:,index2), ...
            dist, -1);
        if plot
            candidatePoint = plot3(xCand(1) , xCand(2), xCand(3), 'go');
        end
        nearIndex = check(xCand, x(:,frontiers{k}.inds), 0.9*dist);
        if (nearIndex == 0)
            newIndex = numPts + 1;
            xCand = NewtonOneStep(xCand, f, grad);
            if plot
                candidatePoint2 = plot3(xCand(1) , xCand(2), xCand(3), 'ro');
            end
        
            x(:,newIndex) = xCand;
            gradX(:,newIndex) = grad(x(:,newIndex));
            faces = [faces; [index1, newIndex, index2]];
            
            frontiers{k}.inds = [frontiers{k}.inds newIndex];
            if plot
                plot3(x(1,[index1 newIndex index2]), ...
                    x(2,[index1 newIndex index2]), ...
                    x(3,[index1 newIndex index2]), 'bo-');
                if gradients
                    quiver3(x(1,newIndex) , x(2,newIndex), x(3,newIndex),...
                        gradX(1,newIndex) , gradX(2,newIndex), gradX(3,newIndex),0);
                end
            end
            numPts = newIndex;
            frontiers{k}.numPts = frontiers{k}.numPts + 1;
        elseif nearIndex == 1 || nearIndex == frontiers{k}.numPts
            frontiers{k}.inds = circshift(frontiers{k}.inds,1,2);
            
        elseif nearIndex == 2
            faces = [faces; [index1, frontiers{k}.inds(2), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                    x(2,[index1 frontiers{k}.inds(2)]), ...
                    x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
            end
            frontiers{k}.inds(1) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            
        elseif nearIndex == frontiers{k}.numPts - 1
            faces = [faces; [index1, frontiers{k}.inds(end - 1), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(end - 1)]), ...
                    x(2,[index1 frontiers{k}.inds(end - 1)]), ...
                    x(3,[index1 frontiers{k}.inds(end - 1)]), 'bo-');
            end
            frontiers{k}.inds(end) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            
        else
            numNewFrontiers = numNewFrontiers + 1;
            faces = [faces; [index1, frontiers{k}.inds(nearIndex), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                    x(2,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                    x(3,[index1 frontiers{k}.inds(nearIndex) index2]), 'bo-');
            end
            newFrontiers{numNewFrontiers}.inds = frontiers{k}.inds(nearIndex:end);
            newFrontiers{numNewFrontiers}.numPts = length(newFrontiers{numNewFrontiers}.inds);
            
            frontiers{k}.inds = frontiers{k}.inds(1:nearIndex);
            frontiers{k}.numPts = length(frontiers{k}.inds);
            
        end
        
        if plot
            delete(candidatePoint);
            delete(candidatePoint2);
            delete(activeEdge);
        end
        
    end
    frontiers = [frontiers newFrontiers];
    numFrontiers = numFrontiers + numNewFrontiers;
    newFrontiers = {};
    numNewFrontiers = 0;
    
    if plot
        for k = 1:numFrontiers
            frontiers{k}.plot = plotFrontier(gca, frontiers{k}, x);
        end
        drawnow
        for k = 1:numFrontiers
            delete(frontiers{k}.plot);
        end
    end
    
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