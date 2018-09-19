function []=main()
%%% Variables and definitions
%% A      = sparse matrix representing the discretization of the
%%          Laplacian --
%% nev    = number of eigenvalues - this is the number of occupied states
%% Domain = struct containing information on the physical domain.
%% Atoms  = struct containing information on the atoms.
%% tol    = tolerance parameter for stopping scf iteration.
%% maxtis = maximum number of SCF iterations allowed.
%% fid    = output file id
%%
%% rho    = final charge density found
%% lam    = eigenvalues computed - their number may be larger than nev
%% W      = set of wave functions.
%%============================================================
%%-------------------- include file --

global CG_prec poldeg diagmeth adaptiveScheme

%%%
%%%  Read atoms types and positions
%%%
%%%  Derive everything else!
%%%
%%%   Units
Ry=13.605698066;

% imports a bunch of common flag settings
RSDFTsettings;

% check that variables hold vavlid values
if (CG_prec~=0 && CG_prec~=1)
    CG_prec=0;
end
if (poldeg<1)
    poldeg=1;
end
if (adaptiveScheme~=0 && adaptiveScheme~=1)
    adaptiveScheme=0;
end
if (diagmeth<0 || diagmeth>3)
    diagmeth=3;
end

% Information about each element
elem    = importdata('elements_new.csv');
N_elements=size(elem.data);
%%%%


disp(' ********************** ')
disp(' DATA INPUT FOR RSDFT')
disp(' *********************  ')
disp('------------------------')

in_data=0;
while (in_data~=1 && in_data~=2)
    in_dataString=input('Input data mode: 1 for manual, 2 for file:  ','s');
    in_data=str2double(in_dataString);
    if (isnan(in_data) || (in_data~=1 && in_data~=2))
        disp('    Please enter 1 or 2');
        in_data=0;
    end    
end

disp('   ')
if (in_data==1)
    at=0;
    while (at<1)% must have atleast one element
        atString=input('Input number of different atomic species: ', 's');
        at=str2double(atString);
        if (isnan(at) || at<1)
           disp('    Please enter a number >=1') 
           at=0;
        end    
    end
        
        
    disp('   ')
    for i=1:at
        correctElementName=0;
        while (correctElementName==0)
            typ = input('Element of species, e.g., Mg, C, O, only first 18 elements supported: ', 's' );
            correctElementName=elementExists(typ,elem);
            if (correctElementName==0)
               disp('    Element not recognized'); 
            end    
        end
        
        
        disp('  Coordinates should be in atomic units ')
        disp('  Example atoms at (0,0,0) should be entered as 0 0 0 on each line ')
        disp('  Terminate with /, i.e., 0 0 0 / for the last entry ')
        errorFreeCoordinates=0;
        while (errorFreeCoordinates==0)
            xyz=[];
            readxyz=[' '];
            while isempty(strfind(readxyz,'/'))
                readxyz=[readxyz,' ',input('  Input coordinates ','s')];
            end

            try
                xyz=reshape(str2num(readxyz(1:strfind(readxyz,'/')-1)),3,[])';
                
                if (size(xyz,1)==0)
                    disp('    Error found in entered coordinates, try again')
                    errorFreeCoordinates=0;
                else
                    errorFreeCoordinates=1; 
                end    
            catch
                disp('    Error found in entered coordinates, try again')
                errorFreeCoordinates=0;
            end
            
        end
        
        Z_chargeString=input(' How many electrons should be added/removed from the system? ','s');
        Z_charge=str2double(Z_chargeString);
        if (isnan(Z_charge))
           Z_charge=0; 
        end
        Z_charge=Z_charge*-1;% changing from number of electrons to change in charge, adding an electron adds a negative charge
        
        n_atom(i)=size(xyz,1);
 
        str=struct('typ',typ,'coord',xyz);
        Atoms(i) = str;
    end
elseif (in_data==2)
    fileName = input('What is the name of the file to load from?: ', 's' );
    
    if (numel(fileName)<=4)
        disp('File name is too short, must include file extension');
        return;%stop execution and return to command window
    end
    
    %check in the saved molecules folder
    fileName=horzcat('SavedMolecules',filesep,fileName);
       
    fileExtension=fileName(end-3:end);
    
    if (strcmp(fileExtension,'.dat'))
    
        [fid,message] = fopen(fileName,'r');
        if fid == -1
            disp([message , '. File must be in SavedMolecules folder'])
            return; %stop execution and return to command window
        end

        at = fscanf(fid,'%g');
        for i=1:at
            typ = fscanf(fid,'%s',1);
            readxyz = fscanf(fid,'%g');
            n_atom(i)=length(readxyz)/3;
            xyz = reshape(readxyz,3,[])';
            str=struct('typ',typ,'coord',xyz);
            Atoms(i) = str;
        end
        fclose(fid);
    
    elseif (strcmp(fileExtension,'.mat'))
        
        if (exist(fileName,'file')==0)
            disp('File not found in directory, SavedMolecules');
            return;%stop execution and return to command window
        end
        
        variableNames=who('-file',fileName);
        
        if (any(strcmp(variableNames,'AtomsInMolecule')))
            load(fileName,'AtomsInMolecule');
        else
            disp('File does not contain the correct variable')
            return;%stop execution and return to command window
        end
        Atoms=AtomsInMoleculeToAtomsConverter(AtomsInMolecule);
        at=numel(Atoms);
    else    
       disp('File extension not recognized')
       return; %stop execution and return to command window
    end
    if(exist('Z_charge')==0) Z_charge=0; end
end


for tempIndex=1:numel(Atoms)
    n_atom(tempIndex)=numel(Atoms(tempIndex).coord)/3;
end

if (sum(n_atom)<1)
   disp('No atoms in molecule, exiting RSDFT');
   return;
end    
%%%
%%%%Defaut Technical Parameters
%%%%
% fd_order = 8;                         %% order of finite difference scheme
% maxits   = 40;                         %% max SCF iterations
% tol      = 1.e-03;                     %% tolerance for SCF iteration.
% Fermi_temp  =  500.0 ;                %%  Smear out fermi level
global fd_order maxits tol Fermi_temp

%%%%%%%%%
%%%%%%%%%  Grid spacing from table (takes smallest h)
%%%%%%%%%
%%%%
% This is the number of species that need to be looked at
N_types = length(Atoms);

zelec=0.;
hmin=100.;
for at_typ=1:N_types
    % This gets the first atom and assigns the atomic symbol to typ
    typ = Atoms(at_typ).typ;
    % This iterates through elem looking for the matching element data

    %%% Loop over number of atoms
    for i = 1: N_elements
        if strcmp(typ,elem.textdata{i})
            index=i;
            Z=elem.data(index,2)*n_atom(at_typ);
            h=elem.data(index,4);
            if h<hmin
                hmin=h;
            end
            zelec=zelec+Z;
        end
    end
end

ztest=zelec-Z_charge;
if (ztest < 0)
    disp('    Problem with charge state. Negative number of electrons');
    return
end


%%%%
%%%%  hmin is the smallest recommend grid for a particular atom
%%%%  nev is the number of eigenvalues to be found
%%%%
h=hmin;
num_states=0.6*zelec;
if num_states < 8;
    num_states=8;
end
nev=num_states;

%%%%
%%% Estimate the radius
%%%
% Find ion-ion core repulsion
rmax=0.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_types = length(Atoms);
for at_typ=1:N_types
    typ = Atoms(at_typ).typ;
    %
    for i=1:length(elem.data)
        if strcmp(typ,elem.textdata{i})
            index=i;
            break
        end
    end
    xyz=Atoms(at_typ).coord;
    %retrieve corresponding information from elements.csv
    rsize=elem.data(index,5);
    natoms = size(xyz,1);


    %%%%
    for at1 = 1: natoms
        %%%%
        %% scan all points for atom most removed from the domain center
        %%%%
        xx=xyz(at1,3);
        yy=xyz(at1,2);
        zz=xyz(at1,1);
        rdis=sqrt(xx^2+yy^2+zz^2);
        rs=rdis+rsize;
        if rs>rmax
            rmax=rs;
        end
    end
end
sph_rad=rmax;


fid = fopen('./rsdft.out', 'w');       %% create a output file

nx = fix(2*sph_rad/h) + 1;              %%  Make sure nx is even
nx=2*round(nx/2+0.01);

sph_rad=0.5*h*(nx-1);                  %%  Re-adjust R
ny = nx;                               %% [-radis, radius] in each direction
nz = nx;
Domain = struct('radius', sph_rad,'nx',nx,'ny',ny,'nz',nz,'h', h);

fprintf(fid, ' Number of states: \t%d\n\n', num_states);
fprintf(fid,  'Atom data:\n -------------\n');
fprintf(fid, ' Total # of atom types is %d\n', at);
atom_count = 0;
Atoms(1).coord;
for atm_typ = 1:at
    xyz = Atoms(atm_typ).coord ;
    fprintf(fid, ' There are %d %s atoms\n', length(xyz(:,1)), Atoms(atm_typ).typ);
    fprintf(fid, ' and their coordinates are:\n\n');
    fprintf(fid, '\tx [a.u.]\t\ty [a.u.]\t\tz [a.u.]\n');
    atom_count = atom_count + length(xyz(:,1));
    for i=1:length(xyz(:,1))
        fprintf(fid, '\t%.6f\t\t%.6f\t\t%.6f\n', xyz(i,1), xyz(i,2), xyz(i,3));
    end
    fprintf(fid, '\n');
end
fprintf(fid,' --------------------------------------------------\n');
fprintf(fid,' Total number of atoms :         %d\n\n', atom_count);
fprintf(fid,' Number of states:               %10d \n', num_states);
fprintf(fid,' h grid spacing :                %10.5f  \n', h);
fprintf(fid,' Hamiltonian size :              %10d  \n', nx*ny*nz);
fprintf(fid,' Sphere Radius :                 %10.5f   \n', sph_rad);
fprintf(fid,' # grid points in each direction %10d  \n', nx);
fprintf(fid,' Polynomial degree used :        %10d  \n', poldeg);
fprintf(fid,' Finite difference order :       %10d  \n', fd_order);
fprintf(fid,' --------------------------------------------------\n');

disp('      ')
disp('      ')
disp('******************')
disp('     OUTPUT       ')
disp('******************')
disp('      ')
disp(' Working.....constructing Laplacian matrix...')
%% construct Laplacian operator
A  = (1/(h*h))*fd3d(nx, ny, nz, fd_order);


n = size(A,1);
Hpot = zeros(n,1);
pot = Hpot;
%%
%%
err = 10 + tol;
its = 0;
%%
%%  Set up ion-core potentials/initial screening
%%  Import screening from the atom, i.e., hartree potential and charge.
%%
%%  Ion-Ion repulsion
[E_nuc0]=nuclear(Domain, Atoms);
%disp('Nuclear repulsion in eV')
%disp(E_nuc0*Ry)
%%
%%--------------------  Local part of the pseudopotential
%%
disp(' Working.....setting up ionic potential...')

[rho0, hpot0, Ppot]  = pseudoDiag(Domain, Atoms);

%%%%  If the charge state is not neutral, renomalize
if (Z_charge~=0)
    rho0=(ztest/zelec)*rho0;
    hpot0=(ztest/zelec)*hpot0;
end  

%%%%%%%
hpsum0=sum(rho0.*hpot0)*Ry;

fprintf(fid,' Initial Hartree energy (eV) = %10.5f  \n', hpsum0) ;
%%-------------------- count # atoms for stats
n_atoms  = 0;
for ii = 1:length(Atoms)
    n_atoms = n_atoms + length(Atoms(ii).coord(:,1));
end
%%
%%--------------------  Non-local part of the pseudopotential
%%
disp(' Working.....setting up nonlocal part of ionic potential...')

[vnl] = pseudoNL(Domain, Atoms) ;
%%
%%--------------------  Find a initial guess for the screening potential
%%--------------------  Use screening from Gaussian density
%%
h=Domain.h;
%%
rhoxc = rho0' ./ h^3;
[XCpot,exc] = exc_nspn(Domain, rhoxc, fid);
xcpot=XCpot';
Nelec = nelectrons(Atoms);

if (Z_charge~=0)
   Nelec=Nelec-Z_charge; 
end 

%%
%%------------------- open output file (wfn.dat)
%%
wfnid = fopen('./wfn.dat', 'wb');
%% indx5=length(vnl)
pot = Ppot + hpot0 + 0.5*xcpot;
%%-------------------- SCF LOOP
%%-------------------- when 'preconditioning' is used fall ilu0
if (CG_prec)
    disp('calling ilu0 ...')
    [L, U] = luinc(A,'0');
    disp(' done.')
    PRE = struct('L',L, 'U',U);
end

%%-------------------- clear persistent variables in mixer.
clear mixer;
%%-------------------- SCF LOOP starts here
fprintf(fid, '\n----------------------------------\n\n');

halfAPlusvnl=0.5*A+vnl;% A and vnl do not change in the loop and thus so can this calculation

% initialize variables used in adaptive scheme
if (adaptiveScheme~=0 && sum(n_atom)<=2)
    degreeAdaptiveModifier=.75;
    mAdaptiveModifier=.95;
else
    degreeAdaptiveModifier=1;
    mAdaptiveModifier=1;
end

% SCF Loop
while (err > tol && its <= maxits)
    its = its+1;

    fprintf(1,'  Working ... SCF iter # %d  ... ',its)
    %%-------------------- redefine Hamiltonian
    B =  halfAPlusvnl + spdiags(pot, 0, n, n);
    %%--------------------diagmeth defined in include
    tic;
    if (diagmeth ==1 || (its == 1 && diagmeth == 0))
        disp('calling lanczos..')

        v = randn(n,1);
        [W, lam] = lanczos(B, nev+15, v, nev+(500*mAdaptiveModifier), 1.e-05);
    elseif (its == 1 && diagmeth == 2)
        disp('calling chsubsp..')

        [W, lam] = chsubsp(poldeg*degreeAdaptiveModifier, nev+15, B) ;
        
    elseif (its == 1 && diagmeth == 3)
        disp('calling first_filt...')
        [W, lam] = first_filt(nev+15, B, poldeg) ;
        
    else
        disp('calling chebsf..')

        [W, lam] = chefsi1(W, lam, poldeg*degreeAdaptiveModifier, nev, B) ;
    end
    %%
    diag_time = toc;
    %%---------------------print results
    fprintf(fid,' \n \n SCF iter # %d  ... \n',its);
    fprintf(fid,' Diagonalization time [sec] :\t%f\n\n', diag_time);
    %%-------------------- get occupation factors and fermi-level
    %%  increase fermi temp if does not converge

    [Fermi_level, occup] = occupations(lam(1:nev), Fermi_temp, Nelec, 0.000001);
    %%
    fprintf(fid, '   State  Eigenvalue [Ry]     Eigenvalue [eV]\n\n');
    for i = 1:nev
        eig = lam(i) * 2*Ry;
        ry = eig / Ry;
        occ=occup(i);
        fprintf(fid, '%5d   %15.10f   %18.10f  %5.2f\n', i, ry, eig, occ);
    end
    %%-------------------- get charge density
    %rho = (W(:,1:nev) .* W(:,1:nev)) *2* occup ; changed...
    %(W(:,1:nev)*W(:,1:nev)) to W(:,1:nev).^2
    rho=(W(:,1:nev).^2)*2*occup;

    hrhs = (4*pi/h^3)*(rho-rho0);

    rho=rho/h^3;

    %% trigger timer
    tic;
    if (CG_prec)
        Hpot = pcg (A, hrhs, Hpot, 200, 1.e-04, PRE,'precLU');
    else
        Hpot = pcg (A, hrhs, Hpot, 200, 1.e-04);
    end
    hart_time = toc;
    fprintf(fid, '\nHartree potential time [sec]: \t%f\n\n', hart_time);

    [XCpot,exc] = exc_nspn(Domain, rho, fid);
    potNew = Ppot+0.5*XCpot+Hpot+hpot0;
    
    errNew = norm(potNew - pot) / norm(potNew);
    
    % where parameters are changed as part of the
    % adaptive scheme
    if (adaptiveScheme==0 || errNew>1 || errNew>2*err)
        degreeAdaptiveModifier=1;
        mAdaptiveModifier=1;
    elseif (errNew>err)
        degreeAdaptiveModifier=min(1.1,degreeAdaptiveModifier+.2);
        mAdaptiveModifier=min(1.1,degreeAdaptiveModifier+.05);
    elseif (3*errNew<err)
        degreeAdaptiveModifier=max(.5,degreeAdaptiveModifier-.1);   
        mAdaptiveModifier=max(.9,degreeAdaptiveModifier-.025);
    end
    err=errNew; 
    
    fprintf(fid,'   ... SCF error = %10.2e  \n', err) ;
    fprintf(1,'   ... SCF error = %10.2e  \n', err) ;
    %-------------------- call mixer
    pot = mixer(pot, potNew-pot);
end %% end of while loop
%%%%%%%%%%%%%%%%%%%%%%%%


if (err>tol)
    disp('          ')
    disp('**************************')
    disp(' !!THE SYSTEM DID NOT CONVERGE!!')
    disp('          ')
    disp(' !!THESE ARE THE VALUES FROM THE LAST ITERATION!!')
    disp('**************************')
    disp('         ')
else    
    disp('          ')
    disp('**************************')
    disp(' CONVERGED SOLUTION!! ')
    disp('**************************')
    disp('         ')
end

fprintf(1, '   State  Eigenvalue [Ry]     Eigenvalue [eV]  Occupation \n');
for i = 1:nev
    eig = lam(i) * 2*Ry;
    ry = eig / Ry;
    occ=2*occup(i);
    fprintf(1, '%5d   %15.4f   %18.3f  %10.2f\n', i, ry, eig, occ);
end
%%
%%-------------------- total energy
%%  Sum over eigenvalues--Put everything in Ryd
Esum= sum(lam(1:nev).*occup(1:nev));
Esum0=4*Esum;
%%
%%--------------------   Hartree potential
%%  Factor of two for double counting--converts to Ryd
Hsum0=sum(rho.*(Hpot+hpot0))*h^3;
%%
%%-------------------- Exchange correlaion sum
%% No factor of two because energy is in Ry
%%
Vxcsum0=sum(rho.*XCpot)*h^3;
Excsum0=exc;
%%--------------------  Total electronic energy
%%
E_elec0=  Esum0-Hsum0+Excsum0-Vxcsum0;
%%   Add in nuclear-nuclear repulsion term
E_total0 = E_elec0 + E_nuc0;
%%
%%--------------------  Convert to eV
Esum=Ry*Esum0;
Hsum=Ry*Hsum0;
Excsum=Ry*Excsum0;
E_nuc=Ry*E_nuc0;
E_total=Ry*E_total0;
%%-------------------- actual printing of info on energies..
fprintf(fid,'  \n\n');
fprintf(fid,' Total Energies \n\n');
fprintf(fid,' Sum of eigenvalues      = %10.5f  eV   = ',Esum) ;
fprintf(fid,'  %10.5f  Ry  \n',Esum/Ry) ;
fprintf(fid,' Hartree energy          = %10.5f  eV   = ',Hsum) ;
fprintf(fid,'  %10.5f  Ry  \n',Hsum/Ry) ;
fprintf(fid,' Exchange-corr. energy   = %10.5f  eV   = ',Excsum) ;
fprintf(fid,'  %10.5f  Ry  \n',Excsum/Ry) ;
fprintf(fid,' Ion-ion repulsion       = %10.5f  eV   = ',E_nuc) ;
fprintf(fid,'  %10.5f  Ry  \n',E_nuc/Ry) ;
fprintf(fid,' Total electronic energy = %10.5f  eV   = ',E_total) ;
fprintf(fid,'  %10.5f  Ry  \n',E_total0) ;
fprintf(fid,' Electronic energy/atom  = %10.5f  eV   = ',E_total/n_atoms) ;
fprintf(fid,'  %10.5f  Ry  \n',E_total0/n_atoms) ;

%%-------------------- display energies..
fprintf(1,'  \n');
fprintf(1,' Total Energies \n\n');
fprintf(1,' Sum of eigenvalues      = %10.5f  eV   = ',Esum) ;
fprintf(1,'  %10.4f  Ry  \n',Esum/Ry) ;
fprintf(1,' Hartree energy          = %10.5f  eV   = ',Hsum) ;
fprintf(1,'  %10.4f  Ry  \n',Hsum/Ry) ;
fprintf(1,' Exchange-corr. energy   = %10.5f  eV   = ',Excsum) ;
fprintf(1,'  %10.4f  Ry  \n',Excsum/Ry) ;
fprintf(1,' Ion-ion repulsion       = %10.5f  eV   = ',E_nuc) ;
fprintf(1,'  %10.4f  Ry  \n',E_nuc/Ry) ;
fprintf(1,' Total electronic energy = %10.5f  eV   = ',E_total) ;
fprintf(1,'  %10.4f  Ry  \n',E_total0) ;
fprintf(1,' Electronic energy/atom  = %10.5f  eV   = ',E_total/n_atoms);
fprintf(1,'  %10.4f  Ry  \n',E_total0/n_atoms) ;


%%

%%-------------------- free memory (persistent variables) in mixer.
clear mixer;
%%------------------- Output results
fprintf(fid, '\n Finished\n');

little_big_test = 26;
fwrite(wfnid, little_big_test, 'uint32');

fwrite(wfnid, Domain.radius, 'double');
fwrite(wfnid, Domain.h, 'double');

pot_length = length(pot(:,1));
fwrite(wfnid, pot_length, 'uint32');
fwrite(wfnid, pot, 'double');

rho_length = length(rho(:,1));
fwrite(wfnid, rho_length, 'uint32');
fwrite(wfnid, rho, 'double');

w_length = length(W(:,1));
fwrite(wfnid, w_length, 'uint32' );
fwrite(wfnid, nev, 'uint32');
for i = 1:nev
    fwrite(wfnid, W(:, i), 'double');
end

fclose(wfnid);


fclose(fid); %% Close output file

