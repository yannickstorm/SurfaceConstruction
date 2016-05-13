function updateAngles = updateEdgeAngles(frontier, x, gradX, updateInds)

numUpdates = length(updateInds);
updateAngles = zeros(1, numUpdates);
for i = 1:numUpdates
    if updateInds(i) == 1
        ind = 1;
        ind1 = frontier.numPts;
        ind2 = 2;
    elseif updateInds(i) == frontier.numPts
        ind = frontier.numPts;
        ind1 = frontier.numPts - 1;
        ind2 = 1;
    else
        ind = updateInds(i);
        ind1 = updateInds(i) - 1;
        ind2 = updateInds(i) + 1;
    end
    updateAngles(i) = edgeAngle(...
        x(:,frontier.inds(ind)),...
        x(:,frontier.inds(ind1)),...
        x(:,frontier.inds(ind2)),...
        gradX(:,frontier.inds(ind)));
end