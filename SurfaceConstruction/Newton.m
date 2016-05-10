function x = Newton(xIn, f, grad)

x = xIn;
err = abs(f(x));
while err > 1e-5
    fOut = f(x);
    gradIn = grad(x);
	x = x - fOut/norm(gradIn) * gradIn;
    err = abs(f(x));
end
    
