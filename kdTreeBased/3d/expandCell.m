function [ cellRef ] = expandCell( cellRef, evalFun, validState)
    if(cellRef{1} == 2 && isempty(cellRef{2}))
        xLim = cellRef{3};
        yLim = cellRef{4};
        zLim = cellRef{7};
        centroid = cellRef{5};
		c1 = [sum([xLim(1), centroid(1)])/2, sum([yLim(1), centroid(2)])/2, sum([zLim(1), centroid(3)])/2];
		c2 = [sum([xLim(1), centroid(1)])/2, sum([centroid(2), yLim(2)])/2, sum([zLim(1), centroid(3)])/2];
		c3 = [sum([centroid(1), xLim(2)])/2, sum([yLim(1), centroid(2)])/2, sum([zLim(1), centroid(3)])/2];
		c4 = [sum([centroid(1), xLim(2)])/2, sum([centroid(2), yLim(2)])/2, sum([zLim(1), centroid(3)])/2];
        c5 = [sum([xLim(1), centroid(1)])/2, sum([yLim(1), centroid(2)])/2, sum([centroid(3), zLim(2)])/2];
		c6 = [sum([xLim(1), centroid(1)])/2, sum([centroid(2), yLim(2)])/2, sum([centroid(3), zLim(2)])/2];
		c7 = [sum([centroid(1), xLim(2)])/2, sum([yLim(1), centroid(2)])/2, sum([centroid(3), zLim(2)])/2];
		c8 = [sum([centroid(1), xLim(2)])/2, sum([centroid(2), yLim(2)])/2, sum([centroid(3), zLim(2)])/2];
        
        cellRef{2} = {  {validState,{},[xLim(1), centroid(1)],[yLim(1), centroid(2)],c1, evalFun(c1'),[zLim(1), centroid(3)]},... 
                        {validState,{},[xLim(1), centroid(1)],[centroid(2), yLim(2)],c2, evalFun(c2'),[zLim(1), centroid(3)]},... 
                        {validState,{},[centroid(1), xLim(2)],[yLim(1), centroid(2)],c3, evalFun(c3'),[zLim(1), centroid(3)]},... 
                        {validState,{},[centroid(1), xLim(2)],[centroid(2), yLim(2)],c4, evalFun(c4'),[zLim(1), centroid(3)]},... 
                        {validState,{},[xLim(1), centroid(1)],[yLim(1), centroid(2)],c5, evalFun(c5'),[centroid(3), zLim(2)]},... 
                        {validState,{},[xLim(1), centroid(1)],[centroid(2), yLim(2)],c6, evalFun(c6'),[centroid(3), zLim(2)]},... 
                        {validState,{},[centroid(1), xLim(2)],[yLim(1), centroid(2)],c7, evalFun(c7'),[centroid(3), zLim(2)]},... 
                        {validState,{},[centroid(1), xLim(2)],[centroid(2), yLim(2)],c8, evalFun(c8'),[centroid(3), zLim(2)]}};
    else
        for i = 1:8
           if(cellRef{2}{i}{1} == 2)
               cellRef{2}{i} = expandCell( cellRef{2}{i}, evalFun, validState);
           end
        end
    end
end

