function out = check(xCand, x, dist)
out = 0;

distMat = pdist2(xCand', x');
distMat(1) = inf;
distMat(end) = inf;

[val, ind] = min(distMat);

if val < dist
    out = ind;
end