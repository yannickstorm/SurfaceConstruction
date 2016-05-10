function out = check(xCand, x, dist)
out = 0;

distMat = pdist2(xCand', x');

[val, ind] = min(distMat);

if val < dist
    out = ind;
end