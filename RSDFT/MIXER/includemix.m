%
%   This m-file defines parameters for mixing read by mixers (indicated by
%   mixer.m) that can be used by RSDFT.
%
%   Types of mixing methods:
%   simplemix.m: simple mixing;
%   msecant1.m: multi-secant methods with Type-I update to minimize the
%               change of Jacobian;
%   msecant2.m: multi-secant methods with Type-II update to minimize the
%               change of inverse Jacobian;
%   msecant3.m: hybrid method, that minimizes the change of Jacobian or
%               inverse Jacobian depending on the secant errors.
%   ...
%   See mixer.m for how to choose the mixer.
%

% disp('includmix...'); pause;

% Indicates Broyden-like (EN_like=0) or EN-like (EN_like=1) update.
EN_like = 0;


% group_size is the size of the groups.
% For Broyden's family, set group_size=1 and EN_like=0; using msecant1.m
% gets Broyden's first update; using msecant2.m gets Broyden's second
% update; using msecant3.m gets the hybrid method proposed by Martinez.
% The (first) EN-like algorithm proposed by Yang is obtained by setting
% group_size=1 and EN_like=1 and using msecant1.m .
% If group_size = 0 or > number of available iterates, then all available
% iterates are taken in one group, resulting in the Anderson's family.
% to get (the original) Anderson mixing use msecant2.m and set EN_like=0.
group_size = 1;


% Mixing parameter.
mix = 0.5;


% A hybrid method (msecant3.m) choices Type-I or Type-II update depending
% on the secant errors at each iteration. For the first group of iterates
% the secant errors are undefined; hence preferred type of update has to
% be set.
preferred_type = 1;


% (Relative) tolerance for ill-conditioned linear systems.
tol = eps;


% If the |f_new| is too large relative to |f_old|, then the quadratic model
% or its approximation is not reliable, so perform restart. More precisely,
% |f_old| < restart_factor * |f_new|, then perform restart. In particular,
% if restart_factor is 0, then no restart.
restart_factor = 0.1;
