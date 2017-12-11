function [r] = trid (a, d, c, r) 
%% function [r] = trid (a, d, c, r) 
%% algorithm for solving tridag. systems
%% 
%%  d1  c1              x1     r1
%%  a2  d2  c2          x2     r2
%%      a3  d3  c3      x3  =  r3
%%          a4  d4  c4  x4     r4       
%%              a5  d5  x5     r5
%%
%%  note: a, d, c are all of length n - 
%%  a1 and cn not used. 
%% solves A x = r, where A = tridiag[a, d, c]
%%  with a(1) = c(n) = 0 
%% 
n = length(d);
for k=2:n
    piv = a(k)/d(k-1);
    d(k) = d(k) - piv*c(k-1);
    r(k) = r(k) - piv*r(k-1);
end
r(n) = r(n) / d(n);
for j=n-1:-1:1
    r(j) = (r(j) - c(j)*r(j+1))/d(j);
end
 
