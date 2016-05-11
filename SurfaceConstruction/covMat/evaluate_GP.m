function [f_plus_star] = evaluate_GP(K_inv_f, x_star, PartMeans, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals)

[vmean,vmeander] = ComputePriorValues(x_star,Object,cx,cy,cz,a,b,c,Rot,prior_type);
mu_star = [vmean; vmeander];

%%
if With_normals
    f_plus_star = zeros(size(x_star,1)+1,size(mu_star,2));
else
    f_plus_star = zeros(1,size(mu_star,2));
end

for i=1:size(x_star,2)

    K_star = ComputeKderX1X2(sigma,gamma,x_star(:,i),PartMeans,With_normals);
    f_plus_star(:,i) = K_star*K_inv_f;
    if(mod(i,100) == 0)
        i
    end
end
%%
if ~With_normals
    f_plus_star_bis = zeros(size(x_star,1)+1,size(mu_star,2));
    f_plus_star_bis(1,:) = f_plus_star;
    f_plus_star = f_plus_star_bis;
end
f_plus_star = f_plus_star + mu_star;
end