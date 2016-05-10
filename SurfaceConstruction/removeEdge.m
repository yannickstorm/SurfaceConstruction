function edges = removeEdge(edges, edge)
error('bla')
for i = 1:size(edges,1)
    if edges(i,:) == edge
        edges(i,:) = [];
        return
    end
end