function [v_mean, v_meander] = ComputePriorValues(X,Object,cx,cy,cz,a,b,c,Rot,prior_type)

[D,N] = size(X);
Loc = repmat(Object.Loc,[1 N]);




switch prior_type
    case 'S'
%             A = diag([1/a^2,1/a^2]) ;
%             mean = @(x,Loc)((x - Loc)'* A * (x - Loc) - 1)*a/2;
%             meander = @(x,Loc)(A *(x - Loc))*a;
        %% Sphere:
        if size(Object.Loc,1)==2
            A = diag([1/a^2,1/a^2]) ;
            mean = @(x,Loc)((x - Loc)'* A * (x - Loc) - 1)*a/2;
            meander = @(x,Loc)(A *(x - Loc))*a;
        elseif size(Object.Loc,1)==3
            A = diag([1/a^2,1/a^2,1/a^2]) ;
            mean = @(x,Loc)((x - Loc)'* Rot'* A * Rot * (x - Loc) - 1)*a/2;
            meander = @(x,Loc)(Rot'* A * Rot *(x - Loc))*a;
        end
%         mean = @(x,Loc)((x - Loc)'* A  * (x - Loc) - 1)*a/2;
%         meander = @(x,Loc) (A  *(x - Loc))*a;
%         mean = @(x,Loc)((x - Loc)'* Rot'* A * Rot * (x - Loc) - 1)*a/2;
%         meander = @(x,Loc)(Rot'* A * Rot *(x - Loc))*a;

    case 'E'
        
        %% Ellipsoid:
        A = diag([1/a^2,1/b^2,1/c^2]) ;
%         mean = @(x,Loc)((x - Loc)'* A  * (x - Loc) - 1)*a/2;
%         meander = @(x,Loc) (A  *(x - Loc))*a;
        mean = @(x,Loc)((x - Loc)'* Rot'* A * Rot * (x - Loc) - 1)*a/2;
        meander = @(x,Loc)(Rot'* A * Rot *(x - Loc))*a;

    case 'C'
        
        %% Cylinder:
        R=1;
        A = diag([1 0 1]) ;
        mean = @(x,Loc)(((x - Loc)'* A * (x - Loc) - R^2)/(2 * R));
        meander = @(x,Loc)(((A)*(x - Loc))/R);

    case 'P'
        %% Plane:
        A = [1 0 0]
        A = A/ norm(A);
        mean = @(x,Loc)((A * (x - Loc)  ));
        meander = @(x,Loc)(repmat(A',[1 size(x,2)]));


end

switch prior_type
    
    case {'C', 'E', 'S'}
   
        v_mean = zeros(1,size(X,2));
        parfor i=1:size(X,2)
            v_mean(i) = mean(X(:,i), Object.Loc);
        end
        v_meander = meander(X, Loc);
    case { 'P'}
        
        v_mean = (mean(X, Loc))';
        v_meander = meander(X, Loc);
end