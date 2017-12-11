  function [x_new, m] = simplemix(x1, f1)
%
%   [new_x, m] = simplemix(x1, f1);
%   Mixing parameter (mix) is defined in includemix.m .
%
%   Input:
%   x1 = latest iterate;
%   f1 = f(x1).
%
%   Output:
%   x_new = new estimate of the solution to f(x)=0;
%   m = number of available iterates (always 1 for simple mixing);
%


%% Initialization.
% Include mixing parameters.
clear includemix;  % In case file was changed.
includemix;

%% Simple mixing.
x_new = x1 + mix*f1;
m = 1;
