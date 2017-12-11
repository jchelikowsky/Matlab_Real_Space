function [bound]=lanczosForChefsi1(B, nev, v, m, tol)
%% function [W, lam] = lan(B, nev, v, m)
%% B   = matrix;
%% nev = number of wanted egenvalues;
%% v   = initial vector
%% m   = number of steps
%% tol = tolerance for stopping
%% can do with (full) reorthogonalization (reorth=1)
%% or no -- reorthogonalization (reorth=0)
% this version of lanczos is called by chefsi1 and only does the
% calculations necessary to find bound
reorth = 0;
%%
global enableMexFilesTest

n = size(B,1);
v = v/norm(v);
v1 = v;
v0 = zeros(n,1);
k = 0 ;
bet = 0;
ll = 0;
%--------------------pre-allocating
if (reorth)
    VV   = zeros(n,m);
end

Tmat = zeros(m+1,m+1);
tr   = zeros(1,m);


%%-------------------- main lanczos loop
while (k < m)
    k   = k+1;

    if (reorth)
        VV(:,k) = v1;
    end

    v   =  B*v1 - bet*v0;
    if (enableMexFilesTest==1)
        v2=BTimesv1MinusbetTimesv0(B,v0,v1,bet,findFirstColumnWithNonZeroElement(B));
        if (any(abs(v-v2)>0.000001))
            exception=struct('message','Mex file descreptency for BTimesv1MinusbetTimesv0','identifier',[],'stack',{mfilename('fullpath') 'filler'});
            logError(exception);
        end
    end

    alp = v1'*v;

    if (enableMexFilesTest==0)
        v   = v-alp*v1 ;
    else
        vTemp=v-alp*v1;
        vTemp2=vMinusalpTimesv1(v,v1,alp);
        v=vTemp;
        if (any(abs(vTemp-vTemp2)>0.000001))
            exception=struct('message','Mex file descreptency for vMinusalpTimesv1','identifier',[],'stack',{mfilename('fullpath') 'filler'});
            logError(exception);
        end
    end

    %%-----------------------------
    %% reorth  -- test for reorth. needed!
    %%-----------------------------
    if (reorth)
        t = VV(:,1:k)'*v ;
        v = v - VV(:,1:k)*t ;
    end
    %%-------------------- normalize and store v1
    bet = norm(v);
    v0 = v1;
    v1 = v/bet;

    %%-------------------- update tmat
    Tmat(k,k)   = alp;
    Tmat(k+1,k) = bet;
    Tmat(k,k+1) = bet;
    NTest  = min(5*nev,m) ;     %% when to   start testing
    %%-------------------- tr, ll,  == for plotting
    if ( ( (k >= NTest) && (mod(k,10) == 0 )) || k == m )

        rr  = eig(Tmat(1:k,1:k));
        rr  = sort(rr) ;      %% sort increasingly
        %%
        %%-------------------- upper bound used by chef_si -- not used otherwise.
        %%    bound = max(abs(rr)) + bet ;
        bound = max(abs(rr)) + bet;   %% *max(abs(X1(:,k)));
        tr1 = sum(rr(1:nev));
        ll = ll+1;
        tr(ll) = tr1;
    end
    %
    % stopping criterion based on sum of eigenvalues.
    % make sure this is is well-converged!
    %
    if (ll>1 && (abs(tr(ll)-tr(ll-1)) < tol*tr(ll-1)))
        break
    end
    % end - big while loop
end
