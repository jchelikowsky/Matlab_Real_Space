   function [z,c,d]=fspline(xi,yi)
%% function [z,c,d]=fspline(xi,yi)
%% z's are scaled to avoid unnecessary
%% operations with factor 6. 
  xi = xi(:)';
  yi = yi(:)';
  np1=length(xi);    n = np1-1; 
%
%%-------------------- get hi's 
%
h  = xi(2:np1)-xi(1:n);
%
%%-------------------- get right-hand-side
%
r  = (yi(2:np1)-yi(1:n)) ./ h(1:n);
r  = r(2:n)-r(1:n-1);
%
%%-------------------- set tridiag matrix
%
L = h(1:n-1); 
U = h(2:n) ;
D = 2*(L+U); 
%
%%-------------------- solve system
%
z = trid(L,D,U,r);
%
%%-------------------- set end values to 0
%
z = [0 , z , 0] ;
%% NOTE rhs = divided by 6 -> zi also 
%
%-------------------- finally set ci's, di's
%
c = yi(1:n)./h(1:n)-z(1:n).*h(1:n);
d = yi(2:np1)./h(1:n)-z(2:np1).*h(1:n);

