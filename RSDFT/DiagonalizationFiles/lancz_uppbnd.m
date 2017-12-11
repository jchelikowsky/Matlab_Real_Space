function  [uppbnd, ritzv] = lancz_uppbnd (n, A, k)
%
% Usage: uppbnd  = lancz_uppbnd (n, A, k)
%
% apply k steps Lanczos to get the upper bound of abs(eig(A)).
%  
% Input: 
%        n  --- dimension
%        k  --- (optional) perform k steps of Lanczos
%               if not provided, k =4 (a relatively small k is enough)
%        A  --- (optional) the matrix (or a script name for MV) 
%  
% Output:
%   upperb  --- estimated upper bound for the eigenvalues
%   ritzv   --- ritz values 

%
% -ykz, august 2009, dallas
%
    
    if (nargin < 3), 
        k = 6; 
    else
        k = min(max(k, 6), 10);    %do not go over 10 steps
    end 

    T = zeros(k);
    v = rand(n,1);     
    v = v/norm(v);

    tol = 2.5e-16;    %before ||f|| reaches eps, convergence should
                      %have happened, so no need to ask for a small ||f|| 
    
    upperb = zeros(3,k);  % save the bounds for each step j=2:k

    f = A*v;
    alpha = v'*f;
    f     = f - alpha * v; 
    T(1,1)= alpha;
    beta = norm(f);

    % compute the bounds for j=1, using rayleight quotient and its corresponding r
    upperb(1,1) = alpha + beta;
    upperb(2,1) = upperb(1,1);
    upperb(3,1) = upperb(1,1);
    
    isbreak = 0;  
    
    for j = 2 : k    %run k steps

        if (beta > tol),
            v0 = v;  v = f/beta;   
            f = A*v;
            f = f - v0*beta;
            alpha = v'*f;
            f  = f - v*alpha;
            T(j,j-1) = beta; T(j-1,j) = beta; T(j,j) = alpha;
        else
            isbreak = 1;
            fprintf(' j = %i, invariant subspace found\n', j);
            break
        end

        beta = norm(f);
        if (isbreak ~=1),
            [X, ritzv] = eig(T(1:j,1:j));
        else
            [X, ritzv] = eig(T(1:j-1,1:j-1));
        end
        ritzv = diag(ritzv);
        
        if (beta < 1e-2), 
            beta = beta * 10;
        end
        upperb(1,j) = ritzv(end) + beta;
        upperb(2,j) = ritzv(end) + abs(X(end, end))*beta;
        upperb(3,j) = ritzv(end) + max(abs(X(end, :)))*beta;
        
    end

    uppbnd = (upperb(1,j)+upperb(2,j))/2;
    uppbnd = (upperb(3,j)+uppbnd)/2;

    