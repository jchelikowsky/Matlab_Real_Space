function [y j_out] = fsplevalIO(z,c,d,xi,x,j_in)
%% [y j_out] = fsplevalB(z,c,d,xi,x,j_in)
%% evaluates free spline function at x.
%% x a scalar (only).
%% z, c, d, as output from fspline..
%% j_in  = input value for interval to try first
%%         this interval is [xi(jin) xi(j_in+1)]
%% j_out =  output interval for next call
%%
%% algorithm tests if x is in [x_[j_in] x_[j_in+1] ]
%% if yes -- then computes the value -- if not
%% then finds the correct interval.
%%------------------------------------------

n = length(xi);
%%-------------------- input j_in is incorrect-- restart:
if (j_in <1 || j_in > n-1)
    j_in = 1;
end
%%-------------------- x is outside interval -- bring it to
%%                     one of the boundary points
if (x < xi(1))
    x = xi(1);
    j_in = 1;
elseif (x > xi(n))
    x = xi(n);
    j_in = n-1;
end
j_out = j_in;
%%-------------------- if point not in input interval
%%                     do Binary search
if (not(xi(j_in)<=x && x<=xi(j_in+1)))
    %[j_out fflag]=binary_search(xi,x);
    %% INLINED Binary Search %%%%%%%%%
    n = length(xi);
    fflag = 0;
    ind_low = 1;
    ind_high = n;
    if (x>=xi(ind_low) && x<=xi(ind_high))
        while (ind_high-ind_low>1)
            ind_middle = floor((ind_high+ind_low)/2);
            val_middle = xi(ind_middle);
            if (x<val_middle)
                ind_high = ind_middle;
            else
                ind_low = ind_middle;
            end
        end
        i_int = ind_low;
    else
        fflag = -1;
        i_int= 0;
        %return
    end
    j_out=i_int;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (fflag ~= 0)
        error(' SPLINE ERROR [ in binary search ] ')
    end
end

t1  = xi(j_out+1) - x;
t2  = x - xi(j_out);
h_j_out = xi(j_out+1) - xi(j_out);
y   = t1*(z(j_out)*t1*t1/h_j_out+c(j_out))+t2*(z(j_out+1)*t2*t2/h_j_out+d(j_out));
