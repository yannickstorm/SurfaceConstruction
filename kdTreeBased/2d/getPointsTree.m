function [ points, col ] = getPointsTree( points, col, cell, onlyNeg)
    if(isempty(cell{2}))
        if(cell{1}~=0)
          if(cell{6}(1) > 0)
              if(~onlyNeg)
                 points = [points; cell{5}];
                 col = [col; +1]; 
              end
          else
            points = [points; cell{5}];
            col = [col;-1];
          end
        end
    else
        for i = 1:4
           [points, col] = getPointsTree(points, col, cell{2}{i}, onlyNeg);
        end

    end

end

