function xOut = NewtonOneStep(xIn, f, grad)

fIn = f(xIn);
gradIn = grad(xIn);
xOut = xIn - fIn/norm(gradIn) * gradIn;

    
