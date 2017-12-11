  function [x_new, m] = msecant3(x1, f1)
%
%   Hybrid mult-secant methods for solving nonlinear equations f(x)=0;
%   At each iteration, Type I or II update of approximate inverse Jacobian
%   depending on which one results in the smaller secant error.
%   [new_x, m] = msecant3(x1, f1);
%   Mixing parameters are defined in includemix.m .
%
%   Input:
%   x1 = latest iterate;
%   f1 = f(x1).
%
%   Output:
%   x_new = new estimate of the solution to f(x)=0;
%   m = number of available iterates;
%


%% Initialization.
% Initially the columns of DX and DF store x_{i+1}-x_i and f_{i+1}-f_i,
% respectively. The E and V in the note also share the storage with DX and
% DF, respectively. N stores the last group N_i of the Type I update in the
% note. The last two groups in DX, second to the last group in DF are
% stored in DX1, DX2, and DF2, respectively.
persistent DX DF N DX1 DX2 DF2  group_size mix tol restart_factor EN_stage preferred_type type
% The following parameters are defined in includemix.m:
% group_size: size of groups of secant equations;
% mix: mixing parameter (typically |mix|<=1);
% tol: tolerance for rank-deficient or ill-condition systems;
% restart_factor: When ||old_f||<restart_factor*||new_f||, then restart;
% EN_like: 1 for EN-like update, 0 for Broyden-type update.
% preferred_type: If there is only one group of secant equations, the
% secant error is undefined. In this case preferred_type indicates which
% type of update to use.


m = size(DX,2);  % Number of previous iterates.
% Note that m here is m-1 in the note.
if m == 0  % No previous iterate available in the first iteration.
    % Include mixing parameters.
    clear includemix;  % In case file was changed.
    includemix;
    % Store the current iterate.
    DX(:,1) = x1;
    DF(:,1) = f1;
    % Set the new estimate.
    % x_new = x1 - f1;  % Corresponds to the fixed point iteration.
    x_new = x1 + mix*f1;  % Simple mixing.
    if EN_like == 0  % Broyden-like update.
        EN_stage = 0;
    else  % EN-like update.
        EN_stage = 2;  % Next EN stage.
    end
    return;
end
if group_size == 0
    sz = m+1;  % Take all available iterates in one group.
else
    sz = group_size;
end


%% When the current iterate is bad, perform restart.
% More precisely, if ||f0|| < restart_factor*||f1||, then restart, where f0
% is the previous function value. In particular, restart_factor = 0 implies
% never restart.
if EN_stage~=1 && m>=2 && norm(DF(:,m),2)<restart_factor*norm(f1,2)
    x1 = DX(:,m);
    f1 = DF(:,m);
    DX = [];
    DF = [];
    N = [];
    DX1 = [];
    DX2 = [];
    DF2 = [];
    DX(:,1) = x1;
    DF(:,1) = f1;
    % Set the new estimate.
    % x_new = x1 - f1;  % Corresponds to the fixed point iteration.
    x_new = x1 + mix*f1;  % Simple mixing.
    m = 0;
    return;
end


%% Compute res and ngroup and update DX, DF and DX1.
res = mod(m+sz-1,sz)+1;  % Size of the last group.
ngroup = (m-res)/sz;
% ngroup does not count the last group; the number of groups is ngroup+1.
if EN_stage ~= 1
    DX(:,m) = x1 - DX(:,m);  % dx_m = x_{m+1} - x_m.
    DF(:,m) = f1 - DF(:,m);  % df_m = f_{m+1} - f_m.
    DX1(:,res) = DX(:,m);
    DX(:,m) = DX(:,m) + mix*DF(:,m);
    for i=1:ngroup
        DX(:,m) = DX(:,m)-DX(:,(i-1)*sz+1:i*sz)*(DF(:,(i-1)*sz+1:i*sz)'*DF(:,m));
    end
    % DX(:,m) is now E(:,m) in the note.
end


%% Compute new x.
x_new = x1 + mix*f1;
for i=1:ngroup
    x_new = x_new - DX(:,(i-1)*sz+1:i*sz)*(DF(:,(i-1)*sz+1:i*sz)'*f1);
    % DX(:,...) is E(:,...) in the note.
    % DF(:,...) is V(:,...) in the note.
end
% Now deal with the last group.
if EN_stage==1 && sz==res
    x_new = x_new - DX(:,ngroup*sz+1:m)*(DF(:,ngroup*sz+1:m)'*f1);
    DX(:,m+1) = x1;
    DF(:,m+1) = f1;
    EN_stage = 2;
    return;
end
if EN_stage ~= 1
    N(:,res) = -mix*DX1(:,res);
    for i=1:ngroup
        N(:,res) = N(:,res)+DF(:,(i-1)*sz+1:i*sz)*(DX(:,(i-1)*sz+1:i*sz)'*DX1(:,res));
    end
end
% Determine Type-I or Type-II update.
if EN_stage~=1 && ngroup>0
    % If EN_stage==1, then type was already computed/determined in the
    % previous (sub-)iteration, and does not need recomputed again.
    lhs = norm(DF(:,m-res+1:m)'*DF2(:,sz-res+1:sz),'fro')*norm(N(:,1:res)'*DF(:,m-res+1:m),'fro');
    rhs = norm(DF(:,m-res+1:m)'*DF(:,m-res+1:m),'fro')*norm(DX1(:,1:res)'*DX2(:,sz-res+1:sz),'fro');
    if lhs > rhs
        type = 1;
    elseif lhs < rhs
        type = 2;
    else
        type = preferred_type;
    end
else
    type = preferred_type;
end
% Now update.
if type == 1
    M = N(:,1:res)'*DF(:,m-res+1:m);
    C = qinv(M,tol);  % Quasi-inverse of M by QR factorization.
    if res == sz
        DX2 = DX1;
        DF2 = DF(:,m-res+1:m);
        DF(:,m-res+1:m) = N(:,1:res)*C';
        % DF(:,m-sz+1:m) is now V(:,m-sz+1:m) in the note.
        x_new = x_new - DX(:,m-res+1:m)*(DF(:,m-res+1:m)'*f1);
    else
        x_new = x_new - DX(:,m-res+1:m)*((C*N(:,1:res)')*f1);
    end
else  % type == 2.
    C2 = qinv(DF(:,m-res+1:m),tol);
    % C2 is the quasi-inverse of DF(:,m-res+1:m) by QR factorization.
    C = C2*C2';
    if res == sz
        DX2 = DX1;
        DF2 = DF(:,m-res+1:m);
        DF(:,m-res+1:m) = DF(:,m-res+1:m)*C;
        % DF(:,m-sz+1:m) is now V(:,m-sz+1:m) in the note.
        x_new = x_new - DX(:,m-res+1:m)*(DF(:,m-res+1:m)'*f1);
    else
        x_new = x_new - DX(:,m-res+1:m)*((DF(:,m-res+1:m)*C)'*f1);
    end
end

if EN_stage ~= 2
    DX(:,m+1) = x1;
    DF(:,m+1) = f1;
end
if EN_stage == 1
    EN_stage = 2;  % Next EN stage.
elseif EN_stage == 2
    EN_stage = 1;  % Next EN stage.
end
