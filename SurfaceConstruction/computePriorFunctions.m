function [meanValue, meanGrad] = computePriorFunctions(Prior)

cs2 = @(r) (cos(norm(r)/2)^2);
si2 = @(r) (sin(norm(r)/2)^2);
cs = @(r) (cos(norm(r)/2));
si = @(r) (sin(norm(r)/2));

R = @(r) ([     (r(1)^2 - r(2)^2 - r(3)^2) * si2(r)  + norm(r)^2 * cs2(r)    , 2 * si(r)* (r(1)*r(2)*si(r) + norm(r) * r(3) * cs(r)) , 2 * si(r)* (r(1)*r(3)*si(r) - norm(r) * r(2) * cs(r)) ;
             2 * si(r)* (r(1)*r(2)*si(r) - norm(r) * r(3) * cs(r)) ,  (r(2)^2 - r(3)^2 - r(1)^2) * si2(r)  + norm(r)^2 * cs2(r)    ,  2 * si(r)* (r(2)*r(3)*si(r) + norm(r) * r(1) * cs(r));
             , 2 * si(r)* (r(1)*r(3)*si(r) + norm(r) * r(2) * cs(r)) , 2 * si(r)* (r(2)*r(3)*si(r) - norm(r) * r(1) * cs(r)),   (r(3)^2 - r(1)^2 - r(2)^2) * si2(r)  + norm(r)^2 * cs2(r)  ]/norm(r)^2);

         
switch Prior.type
    case 'Ci'
            Pos = Prior.pos(1:2);
            A = diag([1/Prior.param(1)^2,1/Prior.param(1)^2]) ;
            meanValue = @(x)((x - Pos')'* A * (x - Pos') - 1)*Prior.param(1)/2;
            meanGrad = @(x)(A *(x - Pos'))*Prior.param(1);
            
    case 'S'
        %% Sphere:
            
            A = diag([1/Prior.param(1)^2, 1/Prior.param(1)^2, 1/Prior.param(1)^2]);
            
            meanValue = @(x)(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos') - 1));
            
            meanGrad = @(x)(Prior.param(1) *  R(Prior.rot)' * A * R(Prior.rot) * (x-Prior.pos'));

        
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
        
    case 'B'
         
        A1 = diag([1/Prior.param(1)^2, 1/Prior.param(1)^2, 0]);
        A2 = diag([1/Prior.param(1)^2, 1/Prior.param(1)^2, 1/Prior.param(1)^2]);
        
        mCyl = @(x)(Prior.param(1)/2 * ((x-Prior.pos')'* R(Prior.rot)' * A1 * R(Prior.rot) * (x-Prior.pos') - 1));
        gCyl = @(x)(Prior.param(1) *  R(Prior.rot)' * A1 * R(Prior.rot) * (x-Prior.pos'));
        
        mSph = @(x,signe)(Prior.param(1)/2 * (((x-Prior.pos')'* R(Prior.rot)' - signe * [0 0 Prior.param(3)]) * A2 * (R(Prior.rot) * (x-Prior.pos') - signe * [0; 0; Prior.param(3)]) - 1));
        gSph = @(x,signe)(Prior.param(1) *  R(Prior.rot)' * A2 * (R(Prior.rot) * (x-Prior.pos') - signe * [0; 0; Prior.param(3)]));
 
        cond = @(x,signe1, signe2) (signe1 * [0 0 1] * R(Prior.rot) * (x-Prior.pos')  >  signe2 * Prior.param(3));
        meanValue = @(x)( cond(x,1,1) * mSph(x,1) +  (cond(x,-1,-1) && cond(x,1,-1)) * mCyl(x) + cond(x,-1,1) * mSph(x,-1) );
        meanGrad  = @(x)( cond(x,1,1) * gSph(x,1) +  (cond(x,-1,-1) && cond(x,1,-1)) * gCyl(x) + cond(x,-1,1) * gSph(x,-1) );
        
    case 'cube'
         
        A = diag(Prior.param);
        meanValue = @(x)(x'.^2 * A * x.^2 - 1)/6*Prior.param(1);
        meanGrad = @(x)(A * x.^3)*Prior.param(1);
       
end
end