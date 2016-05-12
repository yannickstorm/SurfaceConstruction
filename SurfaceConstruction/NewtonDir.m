function x = NewtonDir(xIn, f, grad)

dir = grad(xIn);
dir = dir/norm(dir);
x = xIn;
err = abs(f(x));
while err > 1e-14
    fCurrent = f(x);
    gradCurrent = grad(x);
    gradDir = gradCurrent' * dir;
	x = x - fCurrent/gradDir * dir;
    err = abs(f(x))/norm(gradCurrent);
end
    
