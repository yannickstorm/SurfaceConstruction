function out = check(xCand, x, inds, dist)
out = 0;

doubleInd1 = find(inds == inds(1));
doubleIndEnd = find(inds == inds(end));
x = x(:,inds);
distMat = pdist2(xCand', x');
%exclude current edge points
distMat(1) = inf; 
distMat(end) = inf;
%exclude points that coincide with current edge points
if length(doubleInd1) > 1
    distMat(doubleInd1(2:end)) = inf;
end
if length(doubleIndEnd) > 1
    distMat(doubleIndEnd(1:end - 1)) = inf;
end

[val, ind] = min(distMat);

if val < dist
    out = ind;
end