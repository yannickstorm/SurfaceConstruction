function x = NewtonDir(xIn, f_plus)

f_tot = f_plus(xIn);

dir = f_tot(2:end);
dir = dir/norm(dir);
x = xIn;
err = f_tot(1);
i = 0;
while err > 1e-14
    i = i + 1
    f_tot = f_plus(x);
    fCurrent = f_tot(1);
    gradCurrent = f_tot(2:end);
    gradDir = gradCurrent' * dir;
	x = x - fCurrent/gradDir * dir;
    f_tot = f_plus(x);
    gradCurrent = f_tot(2:end);
    err = abs(f_tot(1))/norm(gradCurrent);
    
%     if(i > 10)
%        disp('asdfasdf');
%        break;
%     end
end
    
