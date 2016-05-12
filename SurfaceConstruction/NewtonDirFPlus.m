function [x, fGradCurrent] = NewtonDirFPlus(xIn, fPlus)

fGradCurrent = fPlus(xIn);
dir = fGradCurrent(2:end);
dir = dir/norm(dir);
x = xIn;
err = abs(fGradCurrent(1));
while err > 1e-14
    gradDir = fGradCurrent(2:end)' * dir;
	x = x - fGradCurrent(1)/gradDir * dir;
    fGradCurrent = fPlus(x);
    err = abs(fGradCurrent(1))/norm(fGradCurrent(2:end));
end