function initPoints = reduceInitPoints(initPoints, dist)

critical = true;
while critical
    distMat = pdist2(initPoints', initPoints') + diag(inf(1,size(initPoints,2)));
    [minVal, minInd] = min(distMat(:));
    if minVal < 4 * dist
        [~, criticalCol] = ind2sub(size(distMat), minInd);
        initPoints(:,criticalCol) = [];
    else
        critical = false;
    end
end