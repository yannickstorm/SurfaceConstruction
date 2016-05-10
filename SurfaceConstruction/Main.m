close all
clear all
r = 1;
A = 1/r^2 * eye(3);
f = @(x)(r/2 * (x' * A * x - 1));
grad = @(x)(r * A * x);

% f = @(x)([1 0 0] * x);
% fd = @(x)([1;0;0]);

    
dist = 0.4;
figure
axis equal
hold on
x(:,1) = [r; 0; 0];
x(:,2) = x(:,1) + [0; dist; 0];
x(:,2) = Newton(x(:,2), f, grad);
x(:,3) = thirdPoint(x(:,1), x(:,2), ...
    grad(x(:,1)), grad(x(:,1)), ...
    dist, 1);
x(:,3) = Newton(x(:,3), f, grad);

plot3(x(1,[1:3 1]) , x(2,[1:3 1]), x(3,[1:3 1]), 'bo-','linewidth',2)
plot3(0,0,0,'go')
edges = [1,2;2,3;3,1];
faces = [1 2 3];
j = 1;
numPts = 3;
pause
while j < 2000
    j
    edges(j,:)
    index1 = edges(j,1);
    index2 = edges(j,2);
    pEdge = plot3(x(1,[index1 index2]), x(2,[index1 index2]), x(3,[index1 index2]),...
        'r-','linewidth',2);
    pStart = plot3(x(1,[index1]), x(2,[index1]), x(3,[index1]),...
        'r.','markersize',40);
    
    xCand = thirdPoint(...
        x(:,index1),x(:,index2), ...
        grad(x(:,index1)), grad(x(:,index2)), ...
        dist, -1);
    pCand = plot3(xCand(1), xCand(2), xCand(3), 'ro');
    nearIndex = check(xCand, x, dist);
    if (nearIndex == 0)
        newIndex = numPts + 1;
        x(:,newIndex) = Newton(xCand, f, grad);
        edges = [edges; [index1, newIndex; newIndex, index2]];
        plot3([x(1,index1), x(1,newIndex), x(1,index2)], ...
            [x(2,index1), x(2,newIndex), x(2,index2)], ...
            [x(3,index1), x(3,newIndex), x(3,index2)], 'bo-');
        numPts = newIndex;
        faces = [faces; [index1, newIndex, index2]];
        
    else
        edges = [edges; [index1, nearIndex]];
        plot3([x(1,index1), x(1,nearIndex)], ...
            [x(2,index1), x(2,nearIndex)], ...
            [x(3,index1), x(3,nearIndex)], 'bo-')
        faces = [faces; [index1, nearIndex, index2]];
        edges = removeEdge(edges, [index2, nearIndex]);
        
    end
    j = j + 1;
    pause
    delete(pCand);
    delete(pEdge);
    delete(pStart);
end
patch('faces',faces,'vertices',x','facecolor','green');
% 
% e(2:3,:) = [Index1, IndexN; Index2, IndexN];
% 
% plot3(x(1,IndexN), x(2,IndexN), x(3,IndexN), 'r.')
