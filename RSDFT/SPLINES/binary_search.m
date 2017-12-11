 function [i_int fflag]=binary_search(xi,x)
%% [ind_low fflag] binary_search(xi,x) returns
%% the interval in which x belongs
%% i.e x is such that xi(i_int)<= x <=xi(i_int+1)
%% if the fflag is 0 else returns the fflag -1
%% which means "value not found"

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
    return   
end
%%--------------------

