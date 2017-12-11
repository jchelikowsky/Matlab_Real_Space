function [W, ritzv] = first_filt(nev, H, polm)
%   function [W, ritzv] = first_filt(nev, H, polm) ;
%
% Input:
%   nev   -- number of ritz values to compute (>= # of occupied states)
%   H     -- The hamiltonian matrix
%   polm  -- polynomial degree  (can be optional)
%
% Output:
%   ritzv -- contains the new ritz values
%   W     -- the ortho-normal basis of a approximate subspace
%
    
% -ykz, august 2009, dallas    
    
%------------------------------------------------------------------
global OPTIMIZATIONLEVEL% get access to global variable

DEBUG = 0;
Energ_tol    = 0.08;
if (nargin<3), 
    polm=10; 
else
    if (polm < 10 || polm > 15), polm=10; end
end
max_iter = max(min(floor(60/polm), 5), 3);
n = size(H,1) ;

%
% call upper bound estimator to compute an upper bound
%
[upperb, ritzv] = lancz_uppbnd(n, H);
lowerb = 0.8*ritzv(1)+0.2*ritzv(end-1);

if (DEBUG==1), 
    fprintf('polm=%i, nev=%i, sizeH=%i, max_iter=%i\n', polm, nev, n, max_iter),
    fprintf('upperb=%e,  lowerb=%e\n', upperb, lowerb), 
end

% filter random vectors  (the size of W here saves at least half the memory
% than lanczos or other diagonalization methods)
W = rand(n, nev);

for it = 1 : max_iter

    W = cheb_filter_slim(H, W, polm, lowerb, upperb);

    [W, G] = qr(W,0);
    Vin = H*W;   %the temp Vin can be further reduced to save memory, not done yet.
    [n,n2] = size(Vin);
    if (n2~=nev), error('wrong number of eigenvalues'), end
    if (OPTIMIZATIONLEVEL~=0)
        G=Rayleighritz(Vin,W,n2);
    else
        % reuse G and overwrite it
        G=zeros(n2,n2);
        for j=1:n2
            for i=1:j
                G(i,j) = Vin(:,i)'*W(:,j);
                G(j,i) = G(i,j);
            end
        end
    end

    [Q, D] = eig(G);
    ritzv = diag(D);
    [ritzv indx]  = sort(ritzv) ; %%% sort increasingly
    lowerb  = ritzv(end) ;
    W = W*Q(:,indx) ;
    if (it == 1), 
        tr0 = sum(ritzv(1:nev)); 
    else
        tr1 = sum(ritzv(1:nev)); 
        %diff = abs(tr1-tr0)
        if (abs(tr1 - tr0) <= Energ_tol*abs(tr0))
            break;
        end
        if (DEBUG==1 && it > 1),
            fprintf('iter=%i, lowerb = %e, ritz_diff=%e\n', it, lowerb, abs(tr1-tr0)),
        end
        tr0 = tr1;
    end

    
end

%--------------------------------------------------------------
function [y] = cheb_filter_slim(H, x, polm, low, high)
%   
%  slim Chebshev iteration, non-scaling normalized version
%
  e = (high - low)/2;
  center= (high+low)/2;
  
  y = H*x;
  y = (-y + center*x)/e;
  
  for i = 2: polm
    ynew = H*y;
    ynew = (- ynew + center*y)* 2/e  - x;
    x = y;
    y = ynew;
  end
  
  