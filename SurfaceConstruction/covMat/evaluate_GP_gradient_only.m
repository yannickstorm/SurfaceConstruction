function [f_plus_star] = evaluate_GP_gradient_only(K_inv_f, x_star, PartMeans, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals)

f_plus_star = evaluate_GP(K_inv_f, x_star, PartMeans, sigma, gamma,Object,cx,cy,cz,a,b,c,Rot,prior_type, With_normals);

f_plus_star = [f_plus_star(2:4:end,:);f_plus_star(3:4:end,:);f_plus_star(4:4:end,:)];

end