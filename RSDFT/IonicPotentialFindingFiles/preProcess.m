function I = preProcess(v)
%% function I = preProcess(v) makes a pre-processing on the vector v
%% The objective is to keep only the relevant points of this vector.
%% So we remove the points for which the function is constant... 
%% IN:
%% v      = column vector we need to pre-process
%% 
%% OUT   : 
%% I      = list of indexes 
%%============================================================ 

I = [1];
eps = 10e-6;

for i = 2:length(v)-1
	if (norm(v(i)-v(i-1))>eps)
		I = [I i];
	end
end

I = [I length(v)];
