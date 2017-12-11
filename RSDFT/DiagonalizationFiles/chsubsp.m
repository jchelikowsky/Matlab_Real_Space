function [W, ritzv] = chsubsp(deg, nev, H)
%   function [W, ritzv] = shsubsp(deg, nev, H) ;
%
%   deg   -- polynomial degree
%   nev   -- number of occupied states
%   H     -- The hamiltonian matrix
% out:
%   ritzv -- contains the new ritz values
%   W     -- the approximate invariant subspace
%------------------------------------------------------------------
global OPTIMIZATIONLEVEL enableMexFilesTest% get access to global variable

Lanc_steps   = max(3*nev,450);
%%%changed Max_out to 30 from 10
Energ_tol    = 0.05;
Max_out_iter = 30;
n = size(H,1) ;
%%
%% call Lanczos with Lanc_steps steps for getting the upper interval bound
%%
[W, ritzv, upperb] = lanczosForChsubsp(H, nev, randn(n,1), Lanc_steps) ;
%%
%%-------------------- use previous eigenvalues for lower bound
%%
tr0     = sum(ritzv);
%%
%%--------------------outer iteration
%%
for it = 1:Max_out_iter
    %%-------------------- use the smallest eigv. estimate for scaling.
    lam1    = min(ritzv);
    %--------------------
    lowerb  = max(ritzv) ;
    if (lowerb > upperb), error('bounds are wrong'); end
    %%
    W= ch_filter(H, W, deg, lam1, lowerb, upperb);

    %
    %-------------------- Rayleigh-ritz projection.
    % orthonormalize the basis, should be replaced with better method
    % for real computations
    %
    [W, R] = qr(W,0);
    %-------------------- compute Hhat = W'*H*W
    Vin = H*W;
    [n,n2] = size(Vin);
    if (OPTIMIZATIONLEVEL~=0)
        G=Rayleighritz(Vin,W,n2);
    else
        %preallocating G
        G=zeros(n2,n2);
        for j=1:n2
            for i=1:j
                G(i,j) = Vin(:,i)'*W(:,j);
                G(j,i) = G(i,j);
            end
        end
        if (enableMexFilesTest==1)
            G2=Rayleighritz(Vin,W,n2);
            if (any(abs(G-G2)>0.000001))
                exception=struct('message','Mex file descreptency for Rayleighritz.c','identifier',[],'stack',{mfilename('fullpath') 'filler'});
                logError(exception);
            end
        end
    end



    %%
    [Q, D] = eig(G);
    ritzv = diag(D);
    [ritzv indx]  = sort(ritzv) ; %%% sort increasingly
    W      = W*Q(:,indx) ;
    tr1    = sum(ritzv(1:nev)) ;
    if (abs(tr1 - tr0) < Energ_tol*abs(tr1))
        break;
    end

    tr0 = tr1;
end

%--------------------

