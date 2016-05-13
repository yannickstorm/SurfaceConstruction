function [x, gradX, frontier] = initialiseFrontier(initPoint, fPlus, dist, inds)

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
    dist, 1, sqrt(3)/2);
[xNew,fGrad] = NewtonOneStepFPlus(xCand, fPlus);
x(:,3) = xNew;
gradX(:,3) = fGrad(2:end);

frontier.inds = inds;
frontier.numPts = 3;
frontier.edgeAngles = [pi pi pi];
