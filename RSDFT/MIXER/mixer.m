 function [x_new, m] = mixer(x1, f1)
%
%   function [x_new, m] = mixer(x1, f1)
%

% Replace 'msecant1' (two occurances) by another mixer name to assign the mixer.

persistent clear_mixer

if size(clear_mixer,1)==0
    % disp('clear mixer...'); pause;
    clear msecant1;
    clear_mixer = 0;
end
[x_new, m] = msecant1(x1, f1);
