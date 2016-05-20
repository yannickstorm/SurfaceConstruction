function [faces, vertices] = computeSurface(locations, surfNormals, ...
    Prior, ...
    meanValue, meanGrad, initPoints, dist, plot)
sigma = Prior.Sigma;
gamma = Prior.Gamma;
noiseVals = Prior.noiseVals;
noiseGrad = Prior.noiseGrad;
fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);

covMatData = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad);
RVector = covMatData\fPlusData;

fPlusGP = @(x)(CovMatStar(sigma, gamma, x, locations) * RVector);

fPlus = @(x)([meanValue(x);meanGrad(x)] + fPlusGP(x));

if plot
    gradients = false;
    
%     figure
	axis equal
    hold on
    
    plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
    quiver3(locations(1,:),locations(2,:),locations(3,:),...
        surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');
    

    set(gca,'view',[-116.0000   -2.8000]);
end

initPoints = reduceInitPoints(initPoints, dist);

numInitPoints = size(initPoints, 2)
frontiers = cell(1, numInitPoints);
x = [];
gradX = [];
faces = [];
for k = 1:numInitPoints
    inds = (k - 1) * 3 + (1:3);
    [xNew, gradXNew, frontierNew] = initialiseFrontier(initPoints(:,k), fPlus, dist,inds);
    frontiers{k} = frontierNew;
    x = [x xNew];
    gradX = [gradX gradXNew];
    faces = [faces; inds];    
    if plot
        %     figure
        plot3(x(1,[inds inds(1)]) , x(2,[inds inds(1)]), x(3,[inds inds(1)]), 'bo-')
        if gradients
            quiver3(x(1,inds) , x(2,inds), x(3,inds),...
                gradX(1,inds) , gradX(2,inds), gradX(3,inds),0);
        end
    end
end


numFrontiers = length(frontiers);
frontierPlots = [];

numXPts = 3 * k;
removeFrontiers = [];

jMax = 10000;
j = 1;
while numFrontiers > 0 && j < jMax
    if mod(j,100) == 0
        j
    end
    for k = 1:numFrontiers
        if frontiers{k}.numPts < 3
            continue
        end
        
        index1 = frontiers{k}.inds(end);
        index2 = frontiers{k}.inds(1);
        
        if plot
            activeEdge = plot3(x(1,[index1 index2]), ...
                x(2,[index1 index2]), ...
                x(3,[index1 index2]), 'm-','linewidth',4);
        end
        
        if frontiers{k}.edgeAngles(end) < pi/3
            % Angle at ending edge node small enough to close gap immediately
            faces = [faces; [index1, frontiers{k}.inds(end - 1), index2]];
            if plot
                plot3(x(1,[index2 frontiers{k}.inds(end - 1)]), ...
                    x(2,[index2 frontiers{k}.inds(end - 1)]), ...
                    x(3,[index2 frontiers{k}.inds(end - 1)]), 'bo-');
            end
            frontiers{k}.inds(end) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            frontiers{k}.edgeAngles(end) = [];
            updateInds = [1, frontiers{k}.numPts];
            frontiers{k}.edgeAngles(updateInds) = updateEdgeAngles(frontiers{k}, x, gradX,...
                updateInds);
        elseif frontiers{k}.edgeAngles(1) < pi/3
            % Angle at beginning edge node small enough to close gap immediately
            faces = [faces; [index1, frontiers{k}.inds(2), index2]];
            if plot
                plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                    x(2,[index1 frontiers{k}.inds(2)]), ...
                    x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
            end
            frontiers{k}.inds(1) = [];
            frontiers{k}.numPts = length(frontiers{k}.inds);
            frontiers{k}.edgeAngles(1) = [];
            updateInds = [1, frontiers{k}.numPts];
            frontiers{k}.edgeAngles(updateInds) = updateEdgeAngles(frontiers{k}, x, gradX,...
                updateInds);
        else
            % Generate new candidate point
            xCand = thirdPoint(...
                x(:,index1),x(:,index2), ...
                gradX(:,index1), gradX(:,index2), ...
                dist, -1, 0.5);
            if plot
                candidatePoint = plot3(xCand(1) , xCand(2), xCand(3), 'go');
            end
            nearIndex = check(xCand, x, frontiers{k}.inds, dist);
            if (nearIndex == 0)
                %%%%%%%%%%%%%%%
                % Check for intersection with other frontiers
                intersectWithOther = false;
                for kOther = 1:numFrontiers
                    if kOther ~= k
                        nearIndex = check(xCand, x, frontiers{kOther}.inds, dist);
                        if (nearIndex ~= 0)
                            intersectWithOther = true;
                            break
                        end
                    end
                end
                if intersectWithOther == false
                    % We add a new point to the grid 
                    xCand = thirdPoint(...
                        x(:,index1),x(:,index2), ...
                        gradX(:,index1), gradX(:,index2), ...
                        dist, -1, sqrt(3)/2);
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
                    % Two frontiers are merged
                    faces = [faces; [index1, frontiers{kOther}.inds(nearIndex), index2]];
                    if plot
                        plot3(x(1,[index1 frontiers{kOther}.inds(nearIndex) index2]), ...
                            x(2,[index1 frontiers{kOther}.inds(nearIndex) index2]), ...
                            x(3,[index1 frontiers{kOther}.inds(nearIndex) index2]), 'bo-');
                    end
                    oldEnd = frontiers{k}.numPts;
                    frontiers{k}.inds = [frontiers{k}.inds frontiers{kOther}.inds(nearIndex:end) ...
                        frontiers{kOther}.inds(1:nearIndex)];
                    frontiers{k}.numPts = length(frontiers{k}.inds);
                    updateInds = [1, oldEnd, oldEnd + 1, frontiers{k}.numPts];
                    frontiers{k}.edgeAngles(updateInds) = ...
                        updateEdgeAngles(frontiers{k}, ...
                        x, gradX, updateInds);
                    frontiers{kOther}.numPts = 0;
                end
            elseif nearIndex == frontiers{k}.numPts - 1
                % Candidate point close enough to previous edge point
                faces = [faces; [index1, frontiers{k}.inds(end - 1), index2]];
                if plot
                    plot3(x(1,[index2 frontiers{k}.inds(end - 1)]), ...
                        x(2,[index2 frontiers{k}.inds(end - 1)]), ...
                        x(3,[index2 frontiers{k}.inds(end - 1)]), 'bo-');
                end
                frontiers{k}.inds(end) = [];
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles(end) = [];
                updateInds = [1, frontiers{k}.numPts];
                frontiers{k}.edgeAngles(updateInds) = updateEdgeAngles(frontiers{k}, x, gradX,...
                    updateInds);
            elseif nearIndex == 2
                % Candidate point close enough to next edge point
                faces = [faces; [index1, frontiers{k}.inds(2), index2]];
                if plot
                    plot3(x(1,[index1 frontiers{k}.inds(2)]), ...
                        x(2,[index1 frontiers{k}.inds(2)]), ...
                        x(3,[index1 frontiers{k}.inds(2)]), 'bo-');
                end
                frontiers{k}.inds(1) = [];
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles(1) = [];
                updateInds = [1, frontiers{k}.numPts];
                frontiers{k}.edgeAngles(updateInds) = updateEdgeAngles(frontiers{k}, x, gradX,...
                    updateInds);
            elseif (nearIndex ~= 0)
                % Candidate point close enough to other point on the
                % surface: Split surface
                faces = [faces; [index1, frontiers{k}.inds(nearIndex), index2]];
                if plot
                    plot3(x(1,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                        x(2,[index1 frontiers{k}.inds(nearIndex) index2]), ...
                        x(3,[index1 frontiers{k}.inds(nearIndex) index2]), 'bo-');
                end
                newFrontier.inds = frontiers{k}.inds(nearIndex:end);
                newFrontier.numPts = length(newFrontier.inds);
                newFrontier.edgeAngles = frontiers{k}.edgeAngles(nearIndex:end);
                updateInds = [1, newFrontier.numPts];
                newFrontier.edgeAngles(updateInds) = ...
                    updateEdgeAngles(newFrontier, ...
                    x, gradX, updateInds);
                frontiers = [frontiers newFrontier];
                numFrontiers = numFrontiers + 1;
                frontiers{k}.inds = frontiers{k}.inds(1:nearIndex);
                frontiers{k}.numPts = length(frontiers{k}.inds);
                frontiers{k}.edgeAngles = frontiers{k}.edgeAngles(1:nearIndex);
                updateInds = [1, frontiers{k}.numPts];
                frontiers{k}.edgeAngles(updateInds) = ...
                    updateEdgeAngles(frontiers{k}, ...
                    x, gradX, updateInds);
            end
            if plot
                delete(candidatePoint);
                delete(candidatePoint2);
            end
            
            
        end
        if plot
            frontierPlots = [frontierPlots plotFrontier(gca, frontiers{k}, x)];
            delete(activeEdge);
        end
        
    end
    
    for k = 1:numFrontiers
        if frontiers{k}.numPts < 3
            removeFrontiers = [removeFrontiers k];
        elseif frontiers{k}.numPts == 3
            removeFrontiers = [removeFrontiers k];
            faces = [faces; fliplr(frontiers{k}.inds)];
        end
        
    end
    frontiers(removeFrontiers) = [];
    numFrontiers = numFrontiers - length(removeFrontiers);
    removeFrontiers = [];
    
    
    %%%%%%%%%%%%%%%
    if plot
        drawnow
        delete(frontierPlots);
    end
    
    j = j + 1;
    
end
vertices = x';