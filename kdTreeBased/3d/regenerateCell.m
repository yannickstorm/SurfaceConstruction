function [cell, root] = regenerateCell( cell, X , evalFun, root)
    if(isempty(cell{2})) %If has not childrens check value
        if(cell{1} == 0)
           cell{1} = 2;
           cell = expandCell(cell, evalFun, 1);
           [ cell, root ] = validatePoints( cell, root, evalFun );
        end
    else
        for i=1:8
           if(isInBounds(X, cell{2}{i}{3}, cell{2}{i}{4}))
                [cell{2}{i}, root] = regenerateCell(cell{2}{i}, X, evalFun, root);
           end
        end
    end
end

