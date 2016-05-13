function angle = edgeAngle(x,x1,x2,gradX)
dVec1 = x1 - x;
dVec2 = x2 - x;
if (cross(dVec1,dVec2)' * gradX) < 0
    angle = acos(dVec1' * dVec2/(norm(dVec1) * norm(dVec2)));
else
    angle = pi;
end