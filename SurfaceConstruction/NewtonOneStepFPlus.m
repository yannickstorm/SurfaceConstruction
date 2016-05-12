function [x, fGradCurrent] = NewtonOneStepFPlus(xIn, fPlus)

fGradCurrent = fPlus(xIn);
x = xIn - fGradCurrent(1)/norm(fGradCurrent(2:end)) * fGradCurrent(2:end);
