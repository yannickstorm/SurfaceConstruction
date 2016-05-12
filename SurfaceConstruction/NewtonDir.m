function x = NewtonDir(xIn, f_plus)

f_tot = f_plus(xIn);

dir = f_tot(2:end);
dir = dir/norm(dir);
x = xIn;
err = f_tot(1);
i = 0;
if err < 1e-14
    return;
end
while err > 1e-14
    i = i + 1;
    f_tot = f_plus(x);
    fCurrent = f_tot(1);
    gradCurrent = f_tot(2:end);
    gradDir = gradCurrent' * dir;
% 	x = x - fCurrent/gradDir * dir;
    x = x - fCurrent/norm(gradCurrent) * gradCurrent;
    f_tot = f_plus(x);
    fCurrent = f_tot(1);
    gradCurrent = f_tot(2:end);
    err = abs(f_tot(1))/norm(gradCurrent);
    sig = sign(fCurrent);
    
%     if(i > 10)
%        disp('asdfasdf');
%        break;
%     end
end

gradDir = gradCurrent' * dir;
% x = x - 2 * fCurrent/gradDir * dir;
x = x - 2 * fCurrent/norm(gradCurrent) * gradCurrent;
f_tot = f_plus(x);
fCurrent = f_tot(1);

sig == sign(fCurrent)    
