function [faces, vertices] = computeSurface(locations, surfNormals, ...
    sigma, gamma, noiseVals, noiseGrad, ...
    meanValue, meanGrad, initPoint, dist, plot)

% For bmw_total
% cx = -0.1062;
% cy = -0.0867;
% cz = -1.6419;

% For bmw_11
cx = -6.1655;
cy = -0.0472;
cz = -3.6693;

a = 1.0763;
b = 2.3691;
c = 1.3254;



cs2 = @(r) (cos(norm(r)/2)^2);
si2 = @(r) (sin(norm(r)/2)^2);
cs = @(r) (cos(norm(r)/2));
si = @(r) (sin(norm(r)/2));

R = @(r) ([     (r(1)^2 - r(2)^2 - r(3)^2) * si2(r)  + norm(r)^2 * cs2(r)    , 2 * si(r)* (r(1)*r(2)*si(r) + norm(r) * r(3) * cs(r)) , 2 * si(r)* (r(1)*r(3)*si(r) - norm(r) * r(2) * cs(r)) ;
             2 * si(r)* (r(1)*r(2)*si(r) - norm(r) * r(3) * cs(r)) ,  (r(2)^2 - r(3)^2 - r(1)^2) * si2(r)  + norm(r)^2 * cs2(r)    ,  2 * si(r)* (r(2)*r(3)*si(r) + norm(r) * r(1) * cs(r));
             , 2 * si(r)* (r(1)*r(3)*si(r) + norm(r) * r(2) * cs(r)) , 2 * si(r)* (r(2)*r(3)*si(r) - norm(r) * r(1) * cs(r)),   (r(3)^2 - r(1)^2 - r(2)^2) * si2(r)  + norm(r)^2 * cs2(r)  ]/norm(r)^2);


Rot = [1 0 0; 0 1 0; 0 0 1];%R([-6.0233         0    0.1086]);
Object.Loc = [cx; cy; cz];
prior_type = 'E';
With_normals = 1;

fPlusData = ComputeFplus(locations,surfNormals,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals)

covMatData = ComputeFullKder(sigma,gamma,locations,noiseVals,noiseGrad,With_normals);
RVector = covMatData\fPlusData;

% f = @(x)evaluate_GP_function_only(RVector, x, locations, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals);
grad = @(x)evaluate_GP_gradient_only(RVector, x, locations, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals);
f_plus = @(x)evaluate_GP(RVector, x, locations, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals);


x(:,1) = initPoint;
x(:,1) = NewtonDir(x(:,1), f_plus);
gradX(:,1) = grad(x(:,1));

gradInit = gradX(:,1);
perpVector = cross(gradInit,rand(3,1));
normPerpVector = perpVector/norm(perpVector);
x(:,2) = x(:,1) + dist * normPerpVector;
x(:,2) = NewtonDir(x(:,2), f_plus);
gradX(:,2) = grad(x(:,2));

x(:,3) = thirdPoint(x(:,1), x(:,2), ...
    grad(x(:,1)), grad(x(:,1)), ...
    dist, 1);
x(:,3) = NewtonDir(x(:,3), f_plus);
gradX(:,3) = grad(x(:,3));

if plot
    figure
    axis equal
    hold on
    plot3(x(1,[1:3 1]) , x(2,[1:3 1]), x(3,[1:3 1]), 'bo-')
    plot3(0,0,0,'go')
end

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',30);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',2,'color','r');
xlabel('x');
ylabel('y');
zlabel('z');

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
        xCand = thirdPoint(...
            x(:,index1),x(:,index2), ...
            gradX(:,index1), gradX(:,index2), ...
            dist, -1);
        if plot
            pCand1 = plot3(xCand(1) , xCand(2), xCand(3), 'go');
            pCand = plot3(xCand(1) , xCand(2), xCand(3), 'ro');
        end
        nearIndex = check(xCand, x(:,frontiers{k}.inds), 0.9*dist);
        if (nearIndex == 0)
            newIndex = numPts + 1;
            xCand = NewtonDir(xCand, f_plus);
            x(:,newIndex) = xCand;
            gradX(:,newIndex) = grad(x(:,newIndex));
            faces = [faces; [index1, newIndex, index2]];
            
            frontiers{k}.inds = [frontiers{k}.inds newIndex];
            if plot
                plot3(x(1,[index1 newIndex index2]), ...
                    x(2,[index1 newIndex index2]), ...
                    x(3,[index1 newIndex index2]), 'bo-');
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
            delete(pCand1);
            delete(pCand);
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