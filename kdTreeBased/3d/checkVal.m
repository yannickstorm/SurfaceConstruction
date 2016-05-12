function  val  = checkVal( cell, X )
    val = [nan, 0];
    if(length(cell{2}) == 0) %If has not childrens check value
        val = [cell{6}(1), cell{1}];
    else
        for i=1:8
           if(isInBounds(X, cell{2}{i}{3}, cell{2}{i}{4}, cell{2}{i}{7}))
                val = checkVal(cell{2}{i}, X);
           end
        end
    end
end

