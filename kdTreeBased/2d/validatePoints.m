function [ cell, root ] = validatePoints( cell, root, evalFun )
    if(cell{1} == 1 && isempty(cell{2})) %If has not childrens validate
        incX = cell{3}(2) - cell{3}(1);
        incY = cell{4}(2) - cell{4}(1);
        xup = cell{5} - [0,incY];
        xbot = cell{5} + [0,incY];
        xleft = cell{5} - [incX,0];
        xright = cell{5} + [incX,0];
        
        if(isInBounds(xup, root{3}, root{4}))
            res = checkVal(root, xup);
            valup = sign(res(1));
        else
            valup = sign(cell{6}(1));
        end
        if(isInBounds(xbot, root{3}, root{4}))
            res = checkVal(root, xbot);
            valbot = sign(res(1));
        else
            valbot = sign(cell{6}(1));
        end
        if(isInBounds(xleft, root{3}, root{4}))
            res = checkVal(root, xleft);
            valleft = sign(res(1));
        else
            valleft = sign(cell{6}(1));
        end
        if(isInBounds(xright, root{3}, root{4}))
            res = checkVal(root, xright);
            valright = sign(res(1));
        else
            valright = sign(cell{6}(1));
        end
        
        if(abs(valup + valbot + valleft + valright + sign(cell{6}(1))) == 5)
           cell{1} = 0;
        else
           cell{1} = 2;
        end
        
%         if(cell{6}(1) < 0)
        if(false) % regeneration disabled
           if(valup == 1)
                if(res(2) == 0)
                   [cell2,root] = regenerateCell(root, xup, evalFun, root); 
                end
            end
            if(valbot == 1)
                if(res(2) == 0)
                   [cell2,root] = regenerateCell(root,xbot, evalFun, root); 
                end
            end
            if(valleft == 1)
                if(res(2) == 0)
                   [cell2,root] = regenerateCell(root,xleft, evalFun, root); 
                end
            end
            if(valright == 1)
                if(res(2) == 0)
                   [cell2,root] = regenerateCell(root,xright, evalFun, root); 
                end
            end 
        end  
    else if(~isempty(cell{2}))%If has childrens go deeper
        for i = 1:4
           if(cell{2}{i}{1} ~= 0)
               [cell{2}{i}, root] = validatePoints( cell{2}{i}, root, evalFun);
           end
        end
    end
end

