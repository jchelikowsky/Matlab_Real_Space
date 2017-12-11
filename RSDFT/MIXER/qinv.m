  function B = qinv(A, tol)
%
%   B = qinv(A [, tol]);
%   Quasi-inverse of A by QR factorization with column pivoting.
%   tol is the tolerance, a small positive scalar.
%   The QR factorization of A with column pivoting is denoted by A*P=Q*R,
%   where A is the given m-by-n matrix (m>=n), P is permutation matrix,
%   Q is an m-by-n matrix of orthogonal columns, R is rank-revealing upper
%   triangular matrix in the form [ R11 R12 ; 0 0 ], with R11 nonsingular.
%   In some sense S = [ R11^{-1} 0 ; 0 0 ] plays the role of the inverse
%   of R. We call B = P*S*Q^T the quasi-inverse of B. If A is square and
%   has full-rank, then B is the inverse of A. Note that x=B*b is a
%   solution to the least square problem of minimizing ||Ax-b||_2. Without
%   tolerance this is the solution obtained by A\b in matlab. For better
%   numerical stability one may use quasi-inverse of A, qinv(A) in matlab.
%


if nargin < 2
    tol = eps;
end
[m,n] = size(A);
if m<n
    B = qinv(A',tol)';
    return;
end

tol = tol*norm(A,inf);
[Q,R,p] = qr(A,0);
for r=1:n
    if abs(R(r,r))<tol
        r = r-1;
        break;
    end
end
% R(1:r,1:r) has fully rank numerically.

S = zeros(r,r);
for j=r:-1:1
    S(j,j) = 1/R(j,j);
    for i=j-1:-1:1
        t = 0;
        for k=i+1:r
            t = t+R(i,k)*S(k,j);
        end
        S(i,j) = -t/R(i,i);
    end
end
% S is the inverse of R(1:r,1:r).

B = zeros(n,m);
B(1:r,:) = S*Q(:,1:r)';

% Permute B back to the desired order.
q=zeros(1,n);
for i=1:n
    q(1,p(1,i)) = i;
end
B = B(q,:);

