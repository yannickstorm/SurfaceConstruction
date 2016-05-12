function drawTree( tree )
points=[];
cols=[];

surface = getPointsTree(points, cols, tree, false);

plot3(surface(:,1), surface(:,2),surface(:,3), 'go');

end

