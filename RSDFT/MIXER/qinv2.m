  function B = qinv2(A, tol)
%
%   B = qinv2(A [, tol]);
%   Adapter to invoke pseudo-inverse pinv.m .
%

if nargin < 2
    B = pinv(A);
else
    B = pinv(A,tol);
end
