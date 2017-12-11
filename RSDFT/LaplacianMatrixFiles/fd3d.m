function A = fd3d(nx,ny,nz,degree)
%----------------------------------------------------------------------- 
% function A = fd3d(nx,ny,nz,degree)
% NOTE nx ny and nz must be > 1 for correct matrix generation.
% 
%-----------------------------------------------------------------------
%   
%This was developed when it was noticed kron generates a full matrix and
%then makes it sparse. the approach taken is to generate the appropriate
%diagonals for the neighbor interactions in the x, y and z directions, and
%then use spdiags to generate the matrix in one call.
%
%the following relation is used for neighbor l to create nz by nz blocks of
%ny by ny of blocks of nx by nx matrices (which is a bad description).
%
%x direction: on diagonal l, 
%have repeating blocks of nx-|l| of the coefficient value and |l| zeros
%
%y direction: on diagonal l*nx
%have repeating blocks of nx*(ny-|l|) of the coefficient and nx*|l| zeros
%
%z direction: on diagonal l*nx*ny
%nx*ny*(nz-|l|) of the coefficient are sufficient
%
%note the generated matrix of diagonals, V_diag, actually pads out the
%diagonals to length nx*ny*nz, and flips ordering (rather then regenerate
%the diagonal) because on 'above neighbors' (l > 0) spdiags starts from one 
%end of the diagonal collumn, and an 'below neighbors' (l < 0) spdiags 
%starts from the other end.
%

if degree <= 2
    
    nxyz = nx*ny*nz;
    X1 = kron( ones(ny*nz,1), [ -1*ones(nx-1,1) ; 0 ]);
    Y1 = kron( ones(nz,1), [ -1*ones(nx*(ny-1),1) ; zeros(nx,1) ] );
    V_diag = [ -1*ones(nxyz,1) flipud(Y1) flipud(X1) 6*ones(nxyz,1) X1 Y1 -1*ones(nxyz,1) ];
    N_diag = [ nx*ny nx 1 0 -1 -1*nx -1*nx*ny ];
    A= spdiags( V_diag, N_diag, nxyz, nxyz);
    
elseif degree <= 4 
    
    nxy = nx*ny;
    nxyz = nxy*nz;
    X1 = kron( ones(ny*nz,1), [ (-4/3 * ones(nx-1,1)) ; 0 ]);
    X2 = kron( ones(ny*nz,1), [ (1/12 * ones(nx-2,1)) ; zeros(min(2,nx),1) ]);
    Y1 = kron( ones(nz,1), [ (-4/3 * ones(nx*(ny-1),1)) ; zeros(nx,1) ] );
    Y2 = kron( ones(nz,1), [ (1/12 * ones(nx*(ny-2),1)) ; zeros(min(2*nx,nxy),1) ] );
    V_diag = [ 1/12*ones(nxyz,1) -4/3*ones(nxyz,1) flipud(Y2) flipud(Y1) flipud(X2) flipud(X1) ...
        15/2*ones(nxyz,1) X1 X2 Y1 Y2 -4/3*ones(nxyz,1) 1/12*ones(nxyz,1) ];
    N_diag = [ 2*nxy nxy 2*nx nx 2 1 0 -1 -2 -1*nx -2*nx -1*nxy -2*nxy ];
    A = spdiags( V_diag, N_diag, nxyz, nxyz);
    
elseif degree <= 6 
    
    nxy = nx*ny;
    nxyz = nxy*nz;
    X1 = kron( ones(ny*nz,1), [ (-3/2 * ones(nx-1,1)) ; zeros(min(1,nx),1) ]);
    X2 = kron( ones(ny*nz,1), [ (3/20 * ones(nx-2,1)) ; zeros(min(2,nx),1) ]);
    X3 = kron( ones(ny*nz,1), [ (-1/90 * ones(nx-3,1)) ; zeros(min(3,nx),1) ]);
    Y1 = kron( ones(nz,1), [ (-3/2 * ones(nx*(ny-1),1)) ; zeros(min(nx,nxy),1) ] );
    Y2 = kron( ones(nz,1), [ (3/20 * ones(nx*(ny-2),1)) ; zeros(min(2*nx,nxy),1) ] );
    Y3 = kron( ones(nz,1), [ (-1/90 * ones(nx*(ny-3),1)) ; zeros(min(3*nx,nxy),1) ] );
    V_diag = [ -1/90*ones(nxyz,1) 3/20*ones(nxyz,1) -3/2*ones(nxyz,1) flipud(Y3) flipud(Y2) flipud(Y1) flipud(X3) flipud(X2) flipud(X1) ...
        49/6*ones(nxyz,1) X1 X2 X3 Y1 Y2 Y3 -3/2*ones(nxyz,1) 3/20*ones(nxyz,1) -1/90*ones(nxyz,1) ];
    N_diag = [ 3*nxy 2*nxy nxy 3*nx 2*nx nx 3 2 1 ...
        0 -1 -2 -3 -1*nx -2*nx -3*nx -1*nxy -2*nxy -3*nxy ];
    A = spdiags( V_diag, N_diag, nxyz, nxyz);
  
else

    nxy = nx*ny;
    nxyz = nxy*nz;
    X1 = kron( ones(ny*nz,1), [ (-8/5 * ones(nx-1,1)) ; zeros(min(1,nx),1) ]);
    X2 = kron( ones(ny*nz,1), [ (1/5 * ones(nx-2,1)) ; zeros(min(2,nx),1) ]);
    X3 = kron( ones(ny*nz,1), [ (-8/315 * ones(nx-3,1)) ; zeros(min(3,nx),1) ]);
    X4 = kron( ones(ny*nz,1), [ (1/560 * ones(nx-4,1)) ; zeros(min(4,nx),1) ]);
    Y1 = kron( ones(nz,1), [ (-8/5 * ones(nx*(ny-1),1)) ; zeros(min(nxy,nx),1) ] );
    Y2 = kron( ones(nz,1), [ (1/5 * ones(nx*(ny-2),1)) ; zeros(min(nxy,2*nx),1) ] );
    Y3 = kron( ones(nz,1), [ (-8/315 * ones(nx*(ny-3),1)) ; zeros(min(nxy,3*nx),1) ] );
    Y4 = kron( ones(nz,1), [ (1/560 * ones(nx*(ny-4),1)) ; zeros(min(nxy,4*nx),1) ] );
    V_diag = [ 1/560*ones(nxyz,1) -8/315*ones(nxyz,1) 1/5*ones(nxyz,1) -8/5*ones(nxyz,1) flipud(Y4) flipud(Y3) flipud(Y2) flipud(Y1) flipud(X4) flipud(X3) flipud(X2) flipud(X1) ...
        205/24*ones(nxyz,1) X1 X2 X3 X4 Y1 Y2 Y3 Y4 -8/5*ones(nxyz,1) 1/5*ones(nxyz,1) -8/315*ones(nxyz,1) 1/560*ones(nxyz,1) ];
    N_diag = [ 4*nxy 3*nxy 2*nxy nxy 4*nx 3*nx 2*nx nx 4 3 2 1 ...
        0 -1 -2 -3 -4 -1*nx -2*nx -3*nx -4*nx -1*nxy -2*nxy -3*nxy -4*nxy ];
    A = spdiags( V_diag, N_diag, nxyz, nxyz);
  
end
