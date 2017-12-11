function [x, its] = pcg (A, rhs, x0, m, tol, PRE, precfun)
%% solves A x = rhs by pcg
%%------------------------------------------
global OPTIMIZATIONLEVEL enableMexFilesTest

if (OPTIMIZATIONLEVEL~=0 && nargin==5)
    [x,its]=pcgMexFile(A,rhs,x0,m,tol);
else
    x = x0;
    r = rhs - A * x;
    if (nargin >5)
        z = feval(precfun,PRE,r);
    else
        z = r;
    end
    p = z ;
    ro1 = z' * r;
    tol1 = tol*tol*ro1;
    %%

    its = 0 ;
    while (its < m && ro1 > tol1)
        its = its+1;
        ro = ro1;
        ap = A * p;
        alp = ro / ( ap'* p ) ;
        x = x + alp * p ;
        r = r - alp * ap;
        if (nargin >5)
            %%
            %% Unpreconditioned case
            %%
            z = feval(precfun,PRE,r);
        else
            z = r;
        end
        %%
        ro1= z'*r;
        bet = ro1 / ro ;
        p = z + bet * p;
    end

    if (enableMexFilesTest==1)
        [x2,its2]=pcgMexFile(A,rhs,x0,m,tol);
        if (any(abs(x-x2)>0.000001) || abs(its-its2)>0.000001)
            exception=struct('message','Mex file descreptency for pcgMexFile','identifier',[],'stack',{mfilename('fullpath') 'filler'});
            logError(exception);
        end
    end
    
end
