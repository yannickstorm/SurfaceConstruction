function xOut = NewtonOneStep(xIn, f_plus)
f_tot = f_plus(xIn);
fIn = f_tot(1);
gradIn = f_tot(2:end);
xOut = xIn - fIn/norm(gradIn) * gradIn;

    
