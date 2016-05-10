function x = NewtonDir(xIn, f, grad)

dir = grad(xIn);
dir = dir/norm(dir);
x = xIn;
err = abs(f(x));
while err > 1e-5
    fCurrent = f(x);
    gradCurrent = grad(x);
    gradDir = gradCurrent' * dir;
	x = x - fCurrent/gradDir * dir;
    err = abs(f(x));
end
    
