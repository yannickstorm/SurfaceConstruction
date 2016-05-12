function [sigma,gamma,noiseVals,noiseGrad,cx,cy,cz,a,b,c] = ...
    FindHyperparamRed(PartMeans,SurfNormals,sigma0,gamma0,noiseVals0,noiseGrad0,cx0,cy0,cz0,a0,b0,c0,Optim_var)


sigma = sigma0;
gamma = gamma0;
noiseVals = noiseVals0;
noiseGrad = noiseGrad0;
cx = cx0;
cy = cy0;
cz = cz0;
a = a0;
b = b0;
c = c0;





if sum(Optim_var == [1 0 0])==3 
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,x(1),x(2),noiseVals0,noiseGrad0,cx0,cy0,cz0,a0,b0,c0);
        x0 = [sigma0 gamma0];
elseif sum(Optim_var == [0 1 0])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,sigma0,gamma0,noiseVals0,noiseGrad0,x(1),x(2),x(3),a0,b0,c0);
        x0 = [cx0,cy0,cz0];
elseif sum(Optim_var == [0 0 1])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,sigma0,gamma0,noiseVals0,noiseGrad0,cx0,cy0,cz0,x(1),x(2),x(3));
        x0 = [a0,b0,c0];
elseif sum(Optim_var == [0 1 1])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,sigma0,gamma0,noiseVals0,noiseGrad0,x(1),x(2),x(3),x(4),x(5),x(6));
        x0 = [cx0,cy0,cz0,a0,b0,c0];
elseif sum(Optim_var == [1 0 1])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,x(1),x(2),noiseVals0,noiseGrad0,cx0,cy0,cz0,x(3),x(4),x(5));
        x0 = [sigma0 gamma0,a0,b0,c0];
elseif sum(Optim_var == [1 1 0])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,x(1),x(2),noiseVals0,noiseGrad0,x(3),x(4),x(5),a0,b0,c0);
        x0 = [sigma0 gamma0 cx0,cy0,cz0];
elseif sum(Optim_var == [1 1 1])==3
        NegDataLL = @(x)-DataLogLikelihood(PartMeans,SurfNormals,x(1),x(2),noiseVals0,noiseGrad0,x(3),x(4),x(5),x(6),x(7),x(8));
        x0 = [sigma0 gamma0 cx0,cy0,cz0, a0, b0, c0];
end

% Optimization function call
x = fminsearch(NegDataLL,x0);


if Optim_var(1) 
    sigma = x(1);
    gamma = x(2);
end
if Optim_var(2) 
    cx = x(1 + Optim_var(1)*2);
    cy = x(2 + Optim_var(1)*2);
    cz = x(3 + Optim_var(1)*2);
end
if Optim_var(3)
    a = x(1 + Optim_var(1)*2 + Optim_var(2)*3);
    b = x(2 + Optim_var(1)*2 + Optim_var(2)*3);
    c = x(3 + Optim_var(1)*2 + Optim_var(2)*3);
end
    



