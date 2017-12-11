function [rho0, hpot0, pot]=pseudoDiag(Domain, Atoms)
%
% Set up initial screening and local ionic pseudopotential
%  
%%%  Initial charge density is from the atom.  
%%%  Atomic Hartree potential used as reference.  It is 
%%%  added and subtracted from the total energy.
%%%  Exchange correlation is determined from a superposition 
%%%  of atomic charge
%%%

% parsec_dat;

%%% Creating the splines:
%atoms_data;
splineData;

%%% Opening file -- Right now, this is not used
%fid=fopen('oxy_s_pot.txt','r'); 
%pot_oxy=fscanf(fid,'%g %g', [2 inf]);
%rr_oxy_s=pot_oxy(1,:);
%pot_oxy_s=pot_oxy(2,:);
%n_oxy=length(rr_oxy_s)
%rc_oxy_s=rr_oxy_s(n_oxy)   %% last (largest) value or r 
global OPTIMIZATIONLEVEL enableMexFilesTest% gain access to global variable

%%% Localizing variables
nx   = Domain.nx;
ny   = Domain.ny;
nz   = Domain.nz;
h    = Domain.h;
rad  = Domain.radius;

%%% Pre-allocating memory
ndim = nx*ny*nz;
pot  = zeros(ndim,1);
rho0 = zeros(ndim,1);
hpot0 =zeros(ndim,1);

%%%%
N_types = length(Atoms);
Z_sum=0.;

%%---------------- read elements.csv -- info on atoms
elem=importdata('elements_new.csv') ;

%-------------------- Atom-type loop
for at_typ=1:N_types
    typ = Atoms(at_typ).typ;

    % find the index in elem of the atom element equal to typ
    index=find(strcmp(typ,elem.textdata(:)));
    
    
    xyz=Atoms(at_typ).coord;
%-- retrieving corresponding information from elements.csv
    Z=elem.data(index,2);
    natoms = size(xyz,1);
    Z_sum=Z_sum+Z*natoms;
                      
%%--Searching for atom's data
    for i=1:length(AtomFuncData)
        if strcmp(typ,AtomFuncData(i).atom)
            index=i;
            break;
        end
    end
    
%%--Localizing variables and initializing arrays
% Find out what columns charge,pot_S and hartree are stored in
% note that we assume that the radius is in column 1
   
    i_charge=strmatch('charge',data_list,'exact');
    i_pot_S=strmatch('pot_S',data_list,'exact');
    i_hartree=strmatch('hartree',data_list,'exact');
    
    
    atom_data=AtomFuncData(index).data;
    x_charg=atom_data(:,1);
    y_charg=atom_data(:,i_charge);
    x_pot_s=atom_data(:,1);
    y_pot_s=atom_data(:,i_pot_S);
    x_vhart=atom_data(:,1);
    y_vhart=atom_data(:,i_hartree);
    
%%-------------------- pre-processing the data
    I = preProcess(y_charg);
    x_charg = x_charg(I);
    y_charg = y_charg(I);
        
    I = preProcess(y_pot_s);
    x_pot_s = x_pot_s(I);
    y_pot_s = y_pot_s(I);
        
    I = preProcess(y_vhart);
    x_vhart = x_vhart(I);
    y_vhart = y_vhart(I);
    
%%--Calculating the splines
    [z_chg,c_chg,d_chg]=fspline(x_charg, y_charg);
    [z_p_s,c_p_s,d_p_s]=fspline(x_pot_s, y_pot_s);
    [z_vht,c_vht,d_vht]=fspline(x_vhart, y_vhart);
    
%-- Scan all points for each atom (not optimal but OK)--copied from
%-- ppot.m with a small change of formula to find potential
    for at=1:natoms
        indx = 0;
        
        %-- atom coordinates are given wrt sphere center.. need to adjust
        %-- x_point - x_A = (x_point - x_center) - (x_A -x_center) 
        k=0:nz-1;
        dz = (k.*h - rad - xyz(at,3)).^2;
        
        for k=1:nz
            j=0:ny-1;
            dy=dz(k)+(j.*h-rad-xyz(at,2)).^2;
            
            for j=1:ny
%%		
                i=0:nx-1;
                r1=sqrt(dy(j)+(i.*h - rad  - xyz(at,1)).^2);
                ppot=zeros(nx,1);
                rrho=zeros(nx,1);
                hpot00=zeros(nx,1);
%%----------------- next point in space				 
                %%--------initialization of intervals j_ch;j_p_s and j_vht
                j_ch  = 1;
                j_p_s = 1;
                j_vht = 1;
                

                if (OPTIMIZATIONLEVEL~=0)
                	% a c version of the inner most for loop
                    [ppot,rrho,hpot00]=PsuedoDiagLoops(z_p_s,c_p_s,d_p_s,z_chg,c_chg,d_chg,z_vht,c_vht,d_vht,x_pot_s,x_charg,x_vhart,r1,j_p_s,j_ch,j_vht,nx);
                else
                	for i=1:2:nx
%%----------------- Evaluating the splines
                    	% loop has been unrolled once, can be done because nx
                    	% is always even
                    	iPlusOne=i+1;
                    	[ppot(i,1) j_p_s]   = fsplevalIO(z_p_s,c_p_s,d_p_s,x_pot_s,r1(i),j_p_s);
                    	[rrho(i,1) j_ch]    = fsplevalIO(z_chg,c_chg,d_chg,x_charg,r1(i),j_ch);
                    	[hpot00(i,1) j_vht] = fsplevalIO(z_vht,c_vht,d_vht,x_vhart,r1(i),j_vht);
                    
                    	[ppot(iPlusOne,1) j_p_s]   = fsplevalIO(z_p_s,c_p_s,d_p_s,x_pot_s,r1(iPlusOne),j_p_s);
                    	[rrho(iPlusOne,1) j_ch]    = fsplevalIO(z_chg,c_chg,d_chg,x_charg,r1(iPlusOne),j_ch);
                    	[hpot00(iPlusOne,1) j_vht] = fsplevalIO(z_vht,c_vht,d_vht,x_vhart,r1(iPlusOne),j_vht);                  
%%----------------- done atom-specific calculations - now compute
%%----------------- potentials, charge.    
                	end  %% end x
                end
                
                if (enableMexFilesTest==1)
                    [ppot2,rrho2,hpot002]=PsuedoDiagLoops(z_p_s,c_p_s,d_p_s,z_chg,c_chg,d_chg,z_vht,c_vht,d_vht,x_pot_s,x_charg,x_vhart,r1,j_p_s,j_ch,j_vht,nx);
                    if (any(abs(ppot2-ppot)>0.000001) || any(abs(rrho2-rrho)>0.000001) || any(abs(hpot002-hpot00)>0.000001))
                        exception=struct('message','Mex file descreptency for PsuedoDiagLoops','identifier',[],'stack',{mfilename('fullpath') 'filler'});
                        logError(exception);
                    end    
                end    
                
                rrho=max(0,rrho.*h^3./(4*pi));
                
                indxPlusnx=indx+nx;
                indx = indx+1;
                
                pot(indx:indxPlusnx)   = pot(indx:indxPlusnx)+ppot(:,1);
                rho0(indx:indxPlusnx)  = rho0(indx:indxPlusnx)+rrho(:,1);
                hpot0(indx:indxPlusnx) = hpot0(indx:indxPlusnx)+hpot00(:,1);
                
                indx=indxPlusnx;
            end      %% end y
        end          %% end z
    end              %% end atom
end                  %% end atom_typ

%%---- Check charge
Rsum0=sum(rho0);
rho0=Z_sum*rho0/Rsum0;

