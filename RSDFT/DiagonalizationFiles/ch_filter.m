function [y] = ch_filter(A, x, deg, lam1, low, high)
% function [y] = ch_filter(A, x, deg, lam1, low, high)
%  --> Apply chebyshev filter to x.
%  A    = matrix
%  x    = vector (s) to be filtered
%  lam1 = estimate of lowest eigenvalue - for scaling
%         purposes only [rough estimate OK]
%  [low, high] = interval to be damped.
%
global OPTIMIZATIONLEVEL enableMexFilesTest

e = (high - low)/2;
c = (high+low)/2;
sigma1 = e/(lam1 - c);
sigma  = sigma1;


if (enableMexFilesTest==1)
    %%-------------------- degree 1 term
    y = (A*x - c*x) .* (sigma1/e);
    y2=chebyshevfilterDegree1(A,x,c,sigma1,e);

    if (any(y-y2)>0.000001)
        exception=struct('message','Mex file descreptency for chebyshevfilterDegree1.c','identifier',[],'stack',{mfilename('fullpath') 'filler'});
        logError(exception);
    end

    twoDividedbysigma1=2/sigma1;
    inverseOfe=1/e;

    %%-------------------- loop to degree
    for i = 2: deg
        sigma_new = 1/(twoDividedbysigma1 - sigma);
        t1 = 2*sigma_new*inverseOfe;
        t2 = sigma*sigma_new;

        ynew = (A*y - c*y)*t1 - t2*x;
        ynew2=chebyshevfilterDegreeN(A,x,y,c,t1,t2);

        if (any(ynew-ynew2)>0.000001)
            exception=struct('message','Mex file descreptency for chebyshevfilterDegreeN.c','identifier',[],'stack',{mfilename('fullpath') 'filler'});
            logError(exception);
        end

        x = y;
        y = ynew;
        sigma = sigma_new;
    end
else
    if (OPTIMIZATIONLEVEL~=0)
        %%-------------------- degree 1 term
        y=chebyshevfilterDegree1(A,x,c,sigma1,e);

        twoDividedbysigma1=2/sigma1;
        inverseOfe=1/e;

        %%-------------------- loop to degree
        for i = 2: deg
            sigma_new = 1/(twoDividedbysigma1 - sigma);
            t1 = 2*sigma_new*inverseOfe;
            t2 = sigma*sigma_new;
            ynew=chebyshevfilterDegreeN(A,x,y,c,t1,t2);

            x = y;
            y = ynew;
            sigma = sigma_new;
        end
    else

        %%-------------------- degree 1 term
        y = (A*x - c*x) .* (sigma1/e);

        twoDividedbysigma1=2/sigma1;
        inverseOfe=1/e;
        %%-------------------- loop to degree
        for i = 2: deg
            sigma_new = 1/(twoDividedbysigma1 - sigma);
            t1 = 2*sigma_new*inverseOfe;
            t2 = sigma*sigma_new;
            ynew = (A*y - c*y)*t1 - t2*x;

            x = y;
            y = ynew;
            sigma = sigma_new;
        end

    end

end
