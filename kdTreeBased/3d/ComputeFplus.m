function Fplus = ComputeFplus(PartMeans,SurfNormals,Object,cx,cy,cz,a,b,c)

[D,~] = size(PartMeans);
MeanGrad = ComputeMeanGrad(PartMeans,Object,cx,cy,cz,a,b,c);

MeanGradData(2:D+1,:) = SurfNormals;
Fplus = MeanGradData(:) - MeanGrad;

function Fplus = ComputeMeanGrad(PartMeans,Object,cx,cy,cz,a,b,c)

[D,N] = size(PartMeans);
Loc = repmat(Object.Loc,[1 N]);

Fplus = zeros(D + 1,N);


% mean = @(x,Loc)(1/(2 * R) * (sum((x - Loc).^2,1) - R^2));
% meander = @(x,Loc)(1/R * (x - Loc));

% load Ellipse_plot/semi_axis

% 
% Rx = @(ang)([1 0 0 ; 0 cos(ang) sin(ang); 0, -sin(ang) cos(ang)]);
% Ry = @(ang)([cos(ang)  0 sin(ang) ; 0 1 0 ; -sin(ang) 0 cos(ang)]);
% Rz = @(ang)([cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1]);

A = diag([1/a^2,1/b^2,1/c^2]) ;

mean = @(x,Loc)((x - Loc)'* A * (x - Loc) - 1)*a/2;
meander = @(x,Loc)((A' + A)*(x - Loc))*a/2;

Fplus(1,:) = diag(mean(PartMeans,Loc));
Fplus(2:D+1,:) = meander(PartMeans,Loc);
Fplus = Fplus(:);




