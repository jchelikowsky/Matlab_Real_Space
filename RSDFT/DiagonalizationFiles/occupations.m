 function [c, occup] = occupations(lam, Temp, Nelec, tol)
% function t = occupations(a, b, maxits, tol) 
% bisection algorithm
its = 0;
maxits = 200;
%%%%%%
a= min(lam)-1 ;
%%%%%%%
[fa, dummy] = FermiDirac(lam, a, Temp, Nelec);
lmax=ceil(Nelec/2)+1 ;
b = lam(lmax)+1; 
[fb, dummy] = FermiDirac(lam, b, Temp, Nelec);
% 
 c = (b+a)/2;
 [fc, occup] = FermiDirac(lam, c, Temp, Nelec);

error=2*sum(occup)-Nelec;

if (fa*fb > tol) 
     c = b;
     occup = ones(1:lmax);
     disp(' In bisect - fa*fb > 0') 
     return
end 
while (abs(error) > tol && its < maxits) 
    its = its+1;
    c = (b+a)/2;
     [fc, occup] = FermiDirac(lam, c, Temp, Nelec);
	 error=2*sum(occup)-Nelec;
    if (fc*fb < 0) 
        a = c;
    else
        b = c;
        fb = fc;
    end
end 

%% disp(' electrons count ')
%% disp(Nelec)
%% disp(occup')
%% disp(error)
%% disp(its)




