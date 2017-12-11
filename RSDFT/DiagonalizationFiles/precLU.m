
 function x = precLU(PRE, rhs) 
%% function rhs = precLU(PRE, rhs) 
%% arms preconditioning operation
%% PRE = struct for preconditioner
%%-------------------------------------------------

L = PRE.L;
U = PRE.U; 
x = U \ (L \ rhs);

