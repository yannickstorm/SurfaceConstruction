function [faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist, plot)

fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);

covMatData = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad);
RVector = covMatData\fPlusData;

fPlusGP = @(x)(CovMatStar(sigma, gamma, x, locations) * RVector);

fPlus = @(x)([meanValue(x);meanGrad(x)] + fPlusGP(x));

[x, gradX, frontier] = initialiseFrontier(initPoint, fPlus, dist);
frontiers{1} = frontier;
faces = [1 2 3];

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


numFrontiers = length(frontiers);

numXPts = 3;
removeFrontiers = [];
newFrontiers = {};
numNewFrontiers = 0;

jMax = 10000;
j = 1;
while numFrontiers > 0 && j < jMax
    
    for k = 1:numFrontiers
        if frontiers{k}.numPts == 0
            continue
        end
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
            frontiers{k}.edgeAngles(end) = updateEdgeAngles(frontiers{k}, x, gradX,...
                frontiers{k}.numPts);
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
            frontiers{k}.edgeAngles(1) = updateEdgeAngles(frontiers{k}, x, gradX, 1);
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
                frontiers{k}.edgeAngles(end) = updateEdgeAngles(frontiers{k}, x, gradX,...
                    frontiers{k}.numPts);
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
                frontiers{k}.edgeAngles(1) = updateEdgeAngles(frontiers{k}, x, gradX, 1);
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
                updateInds = [1, newFrontiers{numNewFrontiers}.numPts];
                newFrontiers{numNewFrontiers}.edgeAngles(updateInds) = ...
                    updateEdgeAngles(newFrontiers{numNewFrontiers}, ...
                    x, gradX, updateInds);
                
                frontiers{k}.inds = frontiers{k}.inds(1:nearIndex);
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles = frontiers{k}.edgeAngles(1:nearIndex);
                updateInds = [1, frontiers{k}.numPts];
                frontiers{k}.edgeAngles(updateInds) = ...
                    updateEdgeAngles(frontiers{k}, ...
                    x, gradX, updateInds);
                
            elseif (nearIndex == 0)
                %%%%%%%%%%%%%%%
                % Check for intersection with other frontiers
                intersectWithOther = false;
                for kOther = 1:numFrontiers
                    if kOther ~= k
                        nearIndex = check(xCand, x(:,frontiers{kOther}.inds), 0.9*dist);
                        if (nearIndex ~= 0)
                            intersectWithOther = true;
                            break
                        end
                    end
                end
                if intersectWithOther == false
                    newIndex = numXPts + 1;
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
                    frontiers{k}.numPts = frontiers{k}.numPts + 1;
                    frontiers{k}.edgeAngles = [frontiers{k}.edgeAngles pi];
                    updateInds = [1, frontiers{k}.numPts - 1];
                    frontiers{k}.edgeAngles(updateInds) = ...
                        updateEdgeAngles(frontiers{k}, ...
                        x, gradX, updateInds);
                    if plot
                        plot3(x(1,[index1 newIndex index2]), ...
                            x(2,[index1 newIndex index2]), ...
                            x(3,[index1 newIndex index2]), 'bo-');
                        if gradients
                            quiver3(x(1,newIndex) , x(2,newIndex), x(3,newIndex),...
                                gradX(1,newIndex) , gradX(2,newIndex), gradX(3,newIndex),0);
                        end
                    end
                    numXPts = newIndex;
                else
                    oldEnd = frontiers{k}.numPts;
                    frontiers{k}.inds = [frontiers{k}.inds frontiers{kOther}.inds(nearIndex:end) ...
                        frontiers{kOther}.inds(nearIndex:end) frontiers{kOther}.inds(1:nearIndex)];
                    frontiers{k}.numPts = length(frontiers{k}.inds);
                    updateInds = [1, oldEnd, oldEnd + 1, frontiers{k}.numPts];
                    frontiers{k}.edgeAngles(updateInds) = ...
                        updateEdgeAngles(frontiers{k}, ...
                        x, gradX, updateInds);
                    frontiers{kOther}.numPts = 0;
                end
                
                    
                
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
    
    
    %%%%%%%%%%%%%%%
    % All done, now plot frontiers
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
    
end
vertices = x';