function [meanValue, meanGrad] = computePriorFunctions(Prior)

cs2 = @(r) (cos(norm(r)/2)^2);
si2 = @(r) (sin(norm(r)/2)^2);
cs = @(r) (cos(norm(r)/2));
si = @(r) (sin(norm(r)/2));

R = @(r) ([     (r(1)^2 - r(2)^2 - r(3)^2) * si2(r)  + norm(r)^2 * cs2(r)    , 2 * si(r)* (r(1)*r(2)*si(r) + norm(r) * r(3) * cs(r)) , 2 * si(r)* (r(1)*r(3)*si(r) - norm(r) * r(2) * cs(r)) ;
             2 * si(r)* (r(1)*r(2)*si(r) - norm(r) * r(3) * cs(r)) ,  (r(2)^2 - r(3)^2 - r(1)^2) * si2(r)  + norm(r)^2 * cs2(r)    ,  2 * si(r)* (r(2)*r(3)*si(r) + norm(r) * r(1) * cs(r));
             , 2 * si(r)* (r(1)*r(3)*si(r) + norm(r) * r(2) * cs(r)) , 2 * si(r)* (r(2)*r(3)*si(r) - norm(r) * r(1) * cs(r)),   (r(3)^2 - r(1)^2 - r(2)^2) * si2(r)  + norm(r)^2 * cs2(r)  ]/norm(r)^2);

         
switch Prior.type
    case 'S'
        %% Sphere:
        if size(Prior.pos,2)==2
            A = diag([1/Prior.param(1)^2,1/Prior.param(1)^2]) ;
            meanValue = @(x)((x - Prior.pos)'* A * (x - Prior.pos) - 1)*Prior.param(1)/2;
            meanGrad = @(x)(A *(x - Prior.pos))*Prior.param(1);
        elseif size(Prior.pos,2)==3
            A = diag([1/Prior.param(1)^2, 1/Prior.param(1)^2, 1/Prior.param(1)^2]);
            
            meanValue = @(x)(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos') - 1));
            
            meanGrad = @(x)(Prior.param(1) *  R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos'));

        end
        
    case 'E'
         
        A = diag([1/Prior.param(1)^2, 1/Prior.param(2)^2, 1/Prior.param(3)^2]);
        meanValue = @(x)(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos') - 1));
        meanGrad = @(x)(Prior.param(1) *  R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos'));
        
    case 'C'
         
        A = diag([1/Prior.param(1)^2, 0, 1/Prior.param(3)^2]);
        meanValue = @(x)(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos') - 1));
        meanGrad = @(x)(Prior.param(1) *  R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos'));
    case 'N'
        %% No prior
        meanValue = @(x) Prior.param(1);%(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos') - 1));
        meanGrad = @(x) [0;0;0];%(Prior.param(1) * A * R(Prior.rot) * (x-Prior.pos'));
        
end
end