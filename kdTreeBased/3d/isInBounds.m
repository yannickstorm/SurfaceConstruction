function isIn = isInBounds( X, xLims, yLims, zLims )

    if X(1) > xLims(1) && X(1) < xLims(2)  
        if X(2) > yLims(1) && X(2) < yLims(2)  
            if X(3) > zLims(1) && X(3) < zLims(2) 
                isIn = true;
            else
                isIn = false;
            end
        else
            isIn = false;
        end
    else
        isIn = false;
    end

end

