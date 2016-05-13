function x3 = thirdPoint(x1, x2, grad1, grad2, dist, dir, height)

perp = cross(x2 - x1, (grad1 + grad2)/2);
perp = perp/norm(perp) * dir;
x3 = (x1 + x2)/2 + perp * height * dist;