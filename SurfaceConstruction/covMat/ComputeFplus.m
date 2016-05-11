function Fplus = ComputeFplus(PartMeans,SurfNormals,Object,cx,cy,cz,a,b,c,Rot,prior_type,With_normals)

[D,~] = size(PartMeans);
D = D * With_normals;
MeanGrad = ComputeMeanGrad(PartMeans,Object,cx,cy,cz,a,b,c,Rot,prior_type,With_normals);

if With_normals
    MeanGradData(2:D+1,:) = SurfNormals;
end
Fplus =  - MeanGrad;
if With_normals
    Fplus = MeanGradData(:) - MeanGrad;
end

function Fplus = ComputeMeanGrad(PartMeans,Object,cx,cy,cz,a,b,c,Rot,prior_type,With_normals)

[D,N] = size(PartMeans);
D = D * With_normals;

Fplus = zeros(D + 1,N);

[v_mean, v_meander] = ComputePriorValues(PartMeans,Object,cx,cy,cz,a,b,c,Rot,prior_type);

Fplus(1,:) = v_mean;
if With_normals
    Fplus(2:D+1,:) = v_meander;
end
Fplus = Fplus(:);




