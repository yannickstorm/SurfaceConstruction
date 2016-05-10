function x = Newton(xIn, f, grad)
error('this Newton version needs to be fixed')
x = xIn;
err = abs(f(x));
while err > 1e-5
    fOut = f(x);
    gradIn = grad(x);
	x = x - fOut/norm(gradIn) * gradIn;
    err = abs(f(x));
end
    
