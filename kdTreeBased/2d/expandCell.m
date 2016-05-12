function [ cellRef ] = expandCell( cellRef, evalFun, validState)
    if(cellRef{1} == 2 && isempty(cellRef{2}))
        xLim = cellRef{3};
        yLim = cellRef{4};
        centroid = cellRef{5};
		c1 = [sum([xLim(1), centroid(1)])/2, sum([yLim(1), centroid(2)])/2];
		c2 = [sum([xLim(1), centroid(1)])/2, sum([centroid(2), yLim(2)])/2];
		c3 = [sum([centroid(1), xLim(2)])/2, sum([yLim(1), centroid(2)])/2];
		c4 = [sum([centroid(1), xLim(2)])/2, sum([centroid(2), yLim(2)])/2];
        
        cellRef{2} = {  {validState,{},[xLim(1), centroid(1)],[yLim(1), centroid(2)],c1, evalFun(c1')},... 
                        {validState,{},[xLim(1), centroid(1)],[centroid(2), yLim(2)],c2, evalFun(c2')},... 
                        {validState,{},[centroid(1), xLim(2)],[yLim(1), centroid(2)],c3, evalFun(c3')},... 
                        {validState,{},[centroid(1), xLim(2)],[centroid(2), yLim(2)],c4, evalFun(c4')}};
    else
        for i = 1:4
           if(cellRef{2}{i}{1} == 2)
               cellRef{2}{i} = expandCell( cellRef{2}{i}, evalFun, validState);
           end
        end
    end
end

