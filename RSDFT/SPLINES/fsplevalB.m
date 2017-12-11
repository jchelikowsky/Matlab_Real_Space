function [y] = fsplevalB(z,c,d,xi,x) 
%% [y] fsplevalB(z,c,d,xi,x) 
%% evaluates free spline function at x. 
%% x a scalar (only). 
%% z, c, d, as output from fspline.. 
%%------------------------------------------
n = length(xi);
h(1:n-1)=xi(2:n)-xi(1:n-1);

%% find interval in which x belongs thanks to binary search
[j flag]=binary_search(xi,x);

%% evaluation of the spline function at x
if (flag == 0)
   	t1  = xi(j+1) - x; 
    	t2  = x - xi(j);
    	y   = t1*(z(j)*t1*t1/h(j)+c(j))+t2*(z(j+1)*t2*t2/h(j)+d(j));
else 
    	disp('error : x doesn t belong to xi'); 
end 

return
