function [fe, occup] = FermiDirac(lam, EF, Temp, Nelec) 
%% evaluates fermi-dirac function for given eigenvalues
%%  spin factor = 2
     kT = Temp*6.33327186e-06;
     spin = 1;
     t = 1 + exp( (lam - EF) ./ kT ) ; 
     occup = spin ./ t;
     fe = sum(occup) - Nelec/2; 
%%%--------------------------------------------------- 
