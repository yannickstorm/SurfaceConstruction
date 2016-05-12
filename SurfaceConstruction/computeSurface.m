function [faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist, plot)

fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);

covMatData = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad);
RVector = covMatData\fPlusData;

fPlusGP = @(x)(CovMatStar(sigma, gamma, x, locations) * RVector);

fPlus = @(x)([meanValue(x);meanGrad(x)] + fPlusGP(x));

[xNew,fGrad] = NewtonOneStepFPlus(initPoint, fPlus);
x(:,1) = xNew;
gradX(:,1) = fGrad(2:end);

randVec = [0.8181;
    0.8175;
    0.7224]; % random vector for consistency in results
perpVector = cross(gradX(:,1),randVec);
normPerpVector = perpVector/norm(perpVector);
xCand = x(:,1) + dist * normPerpVector;
[xNew,fGrad] = NewtonOneStepFPlus(xCand, fPlus);
x(:,2) = xNew;
gradX(:,2) = fGrad(2:end);

xCand = thirdPoint(x(:,1), x(:,2), ...
    gradX(:,1), gradX(:,2), ...
    dist, 1);
[xNew,fGrad] = NewtonOneStepFPlus(xCand, fPlus);
x(:,3) = xNew;
gradX(:,3) = fGrad(2:end);

if plot
    %     figure
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
frontiers{1}.edgeAngles = 3;

numFrontiers = 1;
nMax = 10000;
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
        
        if frontiers{k}.edgeAngles(end) < pi/2
            faces = [faces; [index1, frontiers{k}.inds(end - 1), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(end - 1)]), ...
                    x(2,[index1 frontiers{k}.inds(end - 1)]), ...
                    x(3,[index1 frontiers{k}.inds(end - 1)]), 'bo-');
            end
            frontiers{k}.inds(end) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            frontiers{k}.edgeAngles(end) = [];
            frontiers{k}.edgeAngles(end) = edgeAngle(...
                x(:,frontiers{k}.inds(end)),...
                x(:,frontiers{k}.inds(end - 1)),...
                x(:,frontiers{k}.inds(1)),...
                gradX(:,frontiers{k}.inds(end)));
        elseif frontiers{k}.edgeAngles(1) < pi/2
            faces = [faces; [index1, frontiers{k}.inds(2), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                    x(2,[index1 frontiers{k}.inds(2)]), ...
                    x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
            end
            frontiers{k}.inds(1) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            frontiers{k}.edgeAngles(1) = [];
            frontiers{k}.edgeAngles(1) = edgeAngle(...
                x(:,frontiers{k}.inds(1)),...
                x(:,frontiers{k}.inds(end)),...
                x(:,frontiers{k}.inds(2)),...
                gradX(:,frontiers{k}.inds(1)));
        else
            
            xCand = thirdPoint(...
                x(:,index1),x(:,index2), ...
                gradX(:,index1), gradX(:,index2), ...
                dist, -1);
            if plot
                candidatePoint = plot3(xCand(1) , xCand(2), xCand(3), 'go');
            end
            nearIndex = check(xCand, x(:,frontiers{k}.inds), 0.9*dist);
            if (nearIndex == frontiers{k}.numPts - 1)
                faces = [faces; [index1, frontiers{k}.inds(end - 1), index2]];
                if plot
                    plot3(x(1,[index1 frontiers{k}.inds(end - 1)]), ...
                        x(2,[index1 frontiers{k}.inds(end - 1)]), ...
                        x(3,[index1 frontiers{k}.inds(end - 1)]), 'bo-');
                end
                frontiers{k}.inds(end) = [];
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles(end) = [];
                frontiers{k}.edgeAngles(end) = edgeAngle(...
                    x(:,frontiers{k}.inds(end)),...
                    x(:,frontiers{k}.inds(end - 1)),...
                    x(:,frontiers{k}.inds(1)),...
                    gradX(:,frontiers{k}.inds(end)));
            elseif (nearIndex == 2)
                faces = [faces; [index1, frontiers{k}.inds(2), index2]];
                if plot
                    plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                        x(2,[index1 frontiers{k}.inds(2)]), ...
                        x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
                end
                frontiers{k}.inds(1) = [];
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles(1) = [];
                frontiers{k}.edgeAngles(1) = edgeAngle(...
                    x(:,frontiers{k}.inds(1)),...
                    x(:,frontiers{k}.inds(end)),...
                    x(:,frontiers{k}.inds(2)),...
                    gradX(:,frontiers{k}.inds(1)));
            elseif (nearIndex ~= 0)
                numNewFrontiers = numNewFrontiers + 1;
                faces = [faces; [index1, frontiers{k}.inds(nearIndex), index2]];
                if plot
                    plot3(x(1,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                        x(2,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                        x(3,[index1 frontiers{k}.inds(nearIndex) index2]), 'bo-');
                end
                newFrontiers{numNewFrontiers}.inds = frontiers{k}.inds(nearIndex:end);
                newFrontiers{numNewFrontiers}.numPts = length(newFrontiers{numNewFrontiers}.inds);
                newFrontiers{numNewFrontiers}.edgeAngles = frontiers{k}.edgeAngles(nearIndex:end);
                newFrontiers{numNewFrontiers}.edgeAngles(1) = edgeAngle(...
                    x(:,frontiers{numNewFrontiers}.inds(1)),...
                    x(:,frontiers{numNewFrontiers}.inds(end)),...
                    x(:,frontiers{numNewFrontiers}.inds(2)),...
                    gradX(:,frontiers{numNewFrontiers}.inds(1)));
                newFrontiers{numNewFrontiers}.edgeAngles(end) = edgeAngle(...
                    x(:,frontiers{numNewFrontiers}.inds(end)),...
                    x(:,frontiers{numNewFrontiers}.inds(end - 1)),...
                    x(:,frontiers{numNewFrontiers}.inds(1)),...
                    gradX(:,frontiers{numNewFrontiers}.inds(end)));
                frontiers{k}.inds = frontiers{k}.inds(1:nearIndex);
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles = frontiers{k}.edgeAngles(1:nearIndex);
                frontiers{k}.edgeAngles(1) = edgeAngle(...
                    x(:,frontiers{k}.inds(1)),...
                    x(:,frontiers{k}.inds(end)),...
                    x(:,frontiers{k}.inds(2)),...
                    gradX(:,frontiers{k}.inds(1)));
                frontiers{k}.edgeAngles(end) = edgeAngle(...
                    x(:,frontiers{k}.inds(end)),...
                    x(:,frontiers{k}.inds(end - 1)),...
                    x(:,frontiers{k}.inds(1)),...
                    gradX(:,frontiers{k}.inds(end)));
            else
                newIndex = numPts + 1;
                [xNew,fGrad] = NewtonOneStepFPlus(xCand, fPlus);
                x(:,newIndex) = xNew;
                gradX(:,newIndex) = fGrad(2:end);
                if plot
                    candidatePoint2 = plot3(x(1,newIndex) , ...
                        x(2,newIndex) , ...
                        x(3,newIndex) , 'ro');
                end
                faces = [faces; [index1, newIndex, index2]];
                
                frontiers{k}.inds = [frontiers{k}.inds newIndex];
                frontiers{k}.edgeAngles = [frontiers{k}.edgeAngles pi];
                frontiers{k}.edgeAngles(1) = edgeAngle(...
                    x(:,frontiers{k}.inds(1)),...
                    x(:,frontiers{k}.inds(end)),...
                    x(:,frontiers{k}.inds(2)),...
                    gradX(:,frontiers{k}.inds(1)));
                frontiers{k}.edgeAngles(end - 1) = edgeAngle(...
                    x(:,frontiers{k}.inds(end - 1)),...
                    x(:,frontiers{k}.inds(end - 2)),...
                    x(:,frontiers{k}.inds(end)),...
                    gradX(:,frontiers{k}.inds(end - 1)));
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
                
            end
            if plot
                delete(candidatePoint);
                delete(candidatePoint2);
            end
        end
        if plot
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