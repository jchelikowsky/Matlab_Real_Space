function [rhoGrid]=RunRSDFT(Atoms,handles,Z_charge)
%%% Variables and definitions
%% A      = sparse matrix representing the discretization of the
%%          Laplacean --
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
%% Z_charge should be the change in the charge of the molecule, adding 1
% electron should make Z_charge -1

% returns charge density so that it can be visualized
%%============================================================
%%%
%%%  Read atoms types and positions
%%%
%%%  Derive everything else!
%%%
%%%   Units
Ry=13.605698066;

% get variables set in include.m
global CG_prec poldeg diagmeth adaptiveScheme
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


startTotalTime=cputime();


at=numel(Atoms);%at holds the number of elements in the molecule


% this cell array holds all the output to the listbox
% in the ProgessGUI figure
barHandle=barh(handles.ProgressBar,0,0,1);
set(barHandle,'EraseMode','xor');
axis(handles.ProgressBar,[0 1 0 .1]);
axis(handles.ProgressBar, 'off');


listboxArray={'START'};
updateProgressBar(listboxArray,handles,0,barHandle);


n_atom=zeros(1,numel(Atoms));
for tempIndex=1:numel(Atoms)
    n_atom(tempIndex)=numel(Atoms(tempIndex).coord)/3;
end

if (sum(n_atom)<1)
    listboxArray=cellstr('NO ATOMS IN MOLECULE. Ending Calculations');
    updateProgressBar(listboxArray,handles,0,barHandle);
    rhoGrid=[];
    return;
end
%%%
%%%%Defaut Technical Parameters
%%%%
%fd_order = 8;                         %% order of finite difference scheme
%maxits   = 40;                        %% max SCF iterations
%tol      = 1.e-03;                    %% tolerance for SCF iteration.
%Fermi_temp  =  500.0 ;                %%  Smear out fermi level
global fd_order maxits tol Fermi_temp
%%%%%%%%%
%%%%%%%%%  Grid spacing from table (takes smallest h)
%%%%%%%%%
% Information about each element
elem    = importdata('elements_new.csv');
N_elements=size(elem.data);
%%%%
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
    listboxArray=cellstr('Problem with charge state. Negative number of electrons');
    updateProgressBar(listboxArray,handles,0,barHandle);
    rhoGrid=[];
    return
end



%%%%
%%%%  hmin is the smallest recommend grid for a particular atom
%%%%  nev is the number of eigenvalues to be found
%%%%
h=hmin;
num_states=round(0.7*zelec+0.5);
if num_states < 16
    num_states=16;
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

outputFileName='./rsdft.out';
fid = fopen(outputFileName, 'w');       %% create a output file

nx = fix(2*sph_rad/h) + 1;              %%  Make sure nx is even
nx=2*round(nx/2+0.01);
%h  = 2*sph_rad/(nx-1)
sph_rad=0.5*h*(nx-1);                  %%  Re-adjust R
ny = nx;                               %% [-radis, radius] in each direction
nz = nx;
Domain = struct('radius', sph_rad,'nx',nx,'ny',ny,'nz',nz,'h', h);

fprintf(fid, ' Number of states: \t%d\n\n', num_states);
fprintf(fid,  'Atom data:\n -------------\n');
fprintf(fid, ' Total # of atom types is %d\n', at);
atom_count = 0;
Atoms(1).coord ;
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

endFirstInitialization=cputime();

startLaplacianTime=cputime();
listboxArray(numel(listboxArray)+1)=cellstr(' Working.....constructing Laplacian matrix...');
updateProgressBar(listboxArray,handles,0,barHandle);

%% construct Laplacian operator
A  = (1/(h*h))*fd3d(nx, ny, nz, fd_order);


n = size(A,1);
Hpot = zeros(n,1);
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
endLaplacianTime=cputime();

startIonicPotentialTime=cputime();
listboxArray(numel(listboxArray)+1)=cellstr(' Working.....setting up ionic potential...');
updateProgressBar(listboxArray,handles,.03,barHandle);

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
endIonicPotentialTime=cputime();

startNonlocalPotential=cputime();
listboxArray(numel(listboxArray)+1)=cellstr(' Working.....setting up nonlocal part of ionic potential...');
updateProgressBar(listboxArray,handles,.05,barHandle);

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
pot = Ppot + hpot0 + 0.5*xcpot  ;
%%-------------------- SCF LOOP
%%-------------------- when 'preconditioning' is used fall ilu0

endNonlocalPotential=cputime();

if (CG_prec)
    startCGTime=cputime();

    listboxArray(numel(listboxArray)+1)=cellstr('calling ilu0 ...');
    updateProgressBar(listboxArray,handles,.05,barHandle);

    [L, U] = luinc(A,'0');
    %setup.type='nofill';
    %[L U]=ilu(A,setup);

    listboxArray(numel(listboxArray)+1)=cellstr(' done');
    updateProgressBar(listboxArray,handles,.05,barHandle);

    PRE = struct('L',L, 'U',U);

    endCGTime=cputime();
end

%%-------------------- clear persistent variables in mixer.
clear mixer;
%%-------------------- SCF LOOP starts here
fprintf(fid, '\n----------------------------------\n\n');

halfAPlusvnl=0.5*A+vnl;% A and vnl do not change in the loop and thus this calculation an be moved
% outside the loop
percentComplete=.05;

if (adaptiveScheme~=0 && sum(n_atom)<=2)
    degreeAdaptiveModifier=.75;
    mAdaptiveModifier=.95;
else
    degreeAdaptiveModifier=1;
    mAdaptiveModifier=1;
end


totalFirstFiltTime=0;
totalLanczosTime=0;
totalChsubspTime=0;
totalChefsiTime=0;
totalHartreeTime=0;
totalXCPotentialTime=0;
totalDiagTime=0;

totalSCFTime=cputime();

while (err > tol && its <= maxits)
    its = its+1;

    %%-------------------- redefine Hamiltonian
    B =  halfAPlusvnl + spdiags(pot, 0, n, n);
    %%--------------------diagmeth defined in include

    tic;
    if (diagmeth ==1 || (its ==1 && diagmeth == 0))
        startLanczosTime=cputime();
        listboxArray(numel(listboxArray)+1)=cellstr(horzcat('  Working ... SCF iter # ',num2str(its),'  ... ','calling lanczos..'));
        updateProgressBar(listboxArray,handles,percentComplete,barHandle);

        v = randn(n,1);
        [W, lam] = lanczos(B, nev+15, v, nev+(500*mAdaptiveModifier), 1.e-05,percentComplete,handles,barHandle);
        diagcpu = cputime()-startLanczosTime; 
        totalLanczosTime=totalLanczosTime+(diagcpu);
    elseif (its == 1 && diagmeth == 2)
        startChsubspTime=cputime();
        listboxArray(numel(listboxArray)+1)=cellstr(horzcat('  Working ... SCF iter # ',num2str(its),'  ... ','calling chsubsp..'));
        updateProgressBar(listboxArray,handles,percentComplete,barHandle);

        [W, lam] = chsubsp(poldeg*degreeAdaptiveModifier, nev+15, B);
        diagcpu =  cputime()-startChsubspTime;
        
        totalChsubspTime=totalChsubspTime+(diagcpu);
    elseif (its == 1 && diagmeth == 3)
        %
        % added by ykz, it is possible to merge this step with the latter subpsace 
        % filtering so that there will be no first step diagonalization
        %
        startFirstFiltTime=cputime();
        listboxArray(numel(listboxArray)+1)=cellstr(horzcat('  Working ... SCF iter # ',num2str(its),'  ... ','calling first_filt.....'));
        updateProgressBar(listboxArray,handles,percentComplete,barHandle);               
        [W, lam] = first_filt(nev+10, B, poldeg) ;
        diagcpu = cputime()-startFirstFiltTime;
        totalFirstFiltTime=totalFirstFiltTime+(diagcpu);
        
    else
        startChefsiTime=cputime();
        listboxArray(numel(listboxArray)+1)=cellstr(horzcat('  Working ... SCF iter # ',num2str(its),'  ... ','calling chebsf..'));
        updateProgressBar(listboxArray,handles,percentComplete,barHandle);

        [W, lam] = chefsi1(W, lam, poldeg*degreeAdaptiveModifier, nev, B);
        diagcpu = cputime()-startChefsiTime;
        totalChefsiTime=totalChefsiTime+(diagcpu);
    end
    %%
    totalDiagTime = totalDiagTime+diagcpu;
    diag_time = toc;

        
    %%---------------------print results
    fprintf(fid,' \n \n SCF iter # %d  ... \n',its);
    fprintf(fid,' Diagonalization time [sec]: %f,   [cpu sec]: %f\n\n', diag_time, diagcpu);
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
    %rho = (W(:,1:nev) .* W(:,1:nev)) *2* occup ; changed to...
    %(W(:,1:nev)*W(:,1:nev)) to W(:,1:nev).^2
    rho=(W(:,1:nev).^2)*2*occup;

    hrhs = (4*pi/h^3)*(rho-rho0);

    rho=rho./h^3;
    %% trigger timer
    tic;
    if (CG_prec)
        Hpot = pcg (A, hrhs, Hpot, 200, 1.e-04, PRE,'precLU');
    else
        Hpot = pcg (A, hrhs, Hpot, 200, 1.e-04);
    end
    hart_time = toc;
    totalHartreeTime=totalHartreeTime+hart_time;
    fprintf(fid, '\nHartree potential time [sec]: \t%f\n\n', hart_time);

    startEXCTime=cputime();
    [XCpot,exc] = exc_nspn(Domain, rho, fid);
    potNew = Ppot+0.5*XCpot+Hpot+hpot0;
    errNew = norm(potNew - pot) / norm(potNew);
    totalXCPotentialTime=totalXCPotentialTime+(cputime()-startEXCTime);


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

    percentComplete=min(1.05-log10(err/tol)/4,1);

    fprintf(fid,'   ... SCF error = %10.2e  \n', err) ;
    [tempStr]=sprintf('   ... SCF error = %10.2e  ', err) ;
    listboxArray(numel(listboxArray)+1)=cellstr(tempStr);
    %-------------------- call mixer
    pot = mixer(pot, potNew-pot);
end %% end of while loop
%%%%%%%%%%%%%%%%%%%%%%%%
totalSCFTime=cputime()-totalSCFTime;

startOutputTime=cputime();
if (err > tol)
    listboxArray(numel(listboxArray)+1)=cellstr('          ');
    listboxArray(numel(listboxArray)+1)=cellstr('**************************');
    listboxArray(numel(listboxArray)+1)=cellstr(' !!THE SYSTEM DID NOT CONVERGE!!');
    listboxArray(numel(listboxArray)+1)=cellstr(' ');
    listboxArray(numel(listboxArray)+1)=cellstr(' !!THESE ARE THE VALUES FROM THE LAST ITERATION!!');
    listboxArray(numel(listboxArray)+1)=cellstr('**************************');
    listboxArray(numel(listboxArray)+1)=cellstr('          ');
    listboxArray(numel(listboxArray)+1)=cellstr('   State  Eigenvalue [Ry]     Eigenvalue [eV]  Occupation');
else
    listboxArray(numel(listboxArray)+1)=cellstr('          ');
    listboxArray(numel(listboxArray)+1)=cellstr('**************************');
    listboxArray(numel(listboxArray)+1)=cellstr(' CONVERGED SOLUTION!! ');
    listboxArray(numel(listboxArray)+1)=cellstr('**************************');
    listboxArray(numel(listboxArray)+1)=cellstr('          ');
    listboxArray(numel(listboxArray)+1)=cellstr('   State  Eigenvalue [Ry]     Eigenvalue [eV]  Occupation');
end

for i = 1:nev
    eig = lam(i) * 2*Ry;
    ry = eig / Ry;
    occ=2*occup(i);

    [str]=sprintf('%5d   %15.4f      %18.3f           %10.2f', i, ry, eig, occ);
    listboxArray(numel(listboxArray)+1)=cellstr(str);
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
fprintf(fid,'  %10.5f  Ry  \n\n',E_total0/n_atoms) ;

if (err <= tol),
    fprintf(fid,' Self consistency reached within %i steps, diagmeth=%i\n', its, diagmeth) ;
else
    fprintf(fid,' Failed to reach self-consistency within %i steps\n', its) ; 
end
fprintf(fid,' Totel Diagonalization time [sec]: %f \n', totalDiagTime) ;
fprintf(fid,' Totel SCF time [sec]: %f \n', totalSCFTime) ;

%%-------------------- display energies..
listboxArray(numel(listboxArray)+1)=cellstr('          ');

listboxArray(numel(listboxArray)+1)=cellstr(' Total Energies');
listboxArray(numel(listboxArray)+1)=cellstr('          ');

[tempstr]=sprintf(' Sum of eigenvalues         = %10.5f  eV   =   %10.4f  Ry  ',Esum,Esum/Ry) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);

[tempstr]=sprintf(' Hartree energy                = %10.5f  eV   =   %10.4f  Ry  ',Hsum,Hsum/Ry) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);

[tempstr]=sprintf(' Exchange-corr. energy   = %10.5f  eV   =   %10.4f  Ry  ',Excsum,Excsum/Ry) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);

[tempstr]=sprintf(' Ion-ion repulsion              = %10.5f  eV   =   %10.4f  Ry  ',E_nuc,E_nuc/Ry) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);

[tempstr]=sprintf(' Total electronic energy   = %10.5f  eV   =   %10.4f  Ry  ',E_total,E_total0) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);

[tempstr]=sprintf(' Electronic energy/atom   = %10.5f  eV   =   %10.4f  Ry  ',E_total/n_atoms,E_total0/n_atoms) ;
listboxArray(numel(listboxArray)+1)=cellstr(tempstr);



%%

%%-------------------- free memory (persistent variables) in mixer.
clear mixer;
%%------------------- Output results for charge density
fprintf(fid, '\n Finished\n');

little_big_test = nx;
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


listboxArray(numel(listboxArray)+1)=cellstr('DONE');
updateProgressBar(listboxArray,handles,1,barHandle);


rhoGrid=reshape(rho,nx,ny,nz);

fclose(fid); %% Close output file


endTotalTime=cputime();

global profileDataOutput

if (profileDataOutput==0)


    fprintf(1,'Time for intializations before Laplacian- %f seconds\n',endFirstInitialization-startTotalTime);
    fprintf(1,'Time to make Laplacian Matrix- %f seconds\n',endLaplacianTime-startLaplacianTime);
    fprintf(1,'Time for Ionic Potential- %f seconds \n',endIonicPotentialTime-startIonicPotentialTime);
    fprintf(1,'Time for non-local Ionic Potential- %f seconds \n',endNonlocalPotential-startNonlocalPotential);
    if (CG_prec)
        fprintf(1,'Time for preconditioning- %f seconds\n',endCGTime-startCGTime);
    end

    if (diagmeth==0)
        fprintf(1,'Total Lanczos time- %f seconds\n',totalLanczosTime);
        fprintf(1,'Total Chefsi time- %f seconds\n',totalChefsiTime);
        fprintf(1,'    Average Chefsi time- %f seconds\n',totalChefsiTime/(its-1));
    elseif (diagmeth==1)
        fprintf(1,'Total Lanczos time- %f seconds\n',totalLanczosTime);
        fprintf(1,'    Average Lanczos time- %f seconds\n',totalLanczosTime/its);
    else
        fprintf(1,'Total Chsubsp time- %f seconds\n',totalChsubspTime);
        fprintf(1,'Total Chefsi time- %f seconds\n',totalChefsiTime);
        fprintf(1,'    Average Chefsi time- %f seconds\n',totalChefsiTime/(its-1));
    end

    fprintf(1,'Total hartree potential time- %f seconds\n',totalHartreeTime);
    fprintf(1,'    Average hartree potential time- %f seconds\n',totalHartreeTime/its);
    fprintf(1,'Total XC potential time- %f seconds\n',totalXCPotentialTime);
    fprintf(1,'    Average XC potential time- %f seconds\n',totalXCPotentialTime/its);

    fprintf(1,'Time to print output- %f seconds\n',endTotalTime-startOutputTime);

    fprintf(1,'Total time for calculations- %f seconds\n',endTotalTime-startTotalTime);
    fprintf(1,'Number of iterations- %d \n',its);
elseif (profileDataOutput==1)
    dateAndTime=clock();
    dateString=datestr(dateAndTime,21);

    fileID=fopen(horzcat('Profiled-',dateString,'.txt'),'w');

    fprintf(fileID,'Time for intializations before Laplacian= %f seconds\n',endFirstInitialization-startTotalTime);
    fprintf(fileID,'Time to make Laplacian Matrix- %f seconds\n',endLaplacianTime-startLaplacianTime);
    fprintf(fileID,'Time for Ionic Potential- %f seconds \n',endIonicPotentialTime-startIonicPotentialTime);
    fprintf(fileID,'Time for non-local Ionic Potential- %f seconds \n',endNonlocalPotential-startNonlocalPotential);
    if (CG_prec)
        fprintf(fileID,'Time for preconditioning- %f seconds\n',endCGTime-startCGTime);
    end

    if (diagmeth==0)
        fprintf(fileID,'Total Lanczos time- %f seconds\n',totalLanczosTime);
        fprintf(fileID,'Total Chefsi time- %f seconds\n',totalChefsiTime);
        fprintf(fileID,'    Average Chefsi time- %f seconds\n',totalChefsiTime/(its-1));
    elseif (diagmeth==1)
        fprintf(fileID,'Total Lanczos time- %f seconds\n',totalLanczosTime);
        fprintf(fileID,'    Average Lanczos time- %f seconds\n',totalLanczosTime/its);
    else
        fprintf(fileID,'Total Chsubsp time- %f seconds\n',totalChsubspTime);
        fprintf(fileID,'Total Chefsi time- %f seconds\n',totalChefsiTime);
        fprintf(fileID,'    Average Chefsi time- %f seconds\n',totalChefsiTime/(its-1));
    end

    fprintf(fileID,'Total hartree potential time- %f seconds\n',totalHartreeTime);
    fprintf(fileID,'    Average hartree potential time- %f seconds\n',totalHartreeTime/its);
    fprintf(fileID,'Total XC potential time- %f seconds\n',totalXCPotentialTime);
    fprintf(fileID,'    Average XC potential time- %f seconds\n',totalXCPotentialTime/its);

    fprintf(fileID,'Time to print output- %f seconds\n',endTotalTime-startOutputTime);

    fprintf(fileID,'Total time for calculations- %f seconds\n',endTotalTime-startTotalTime);
    fprintf(fileID,'Number of iterations- %d \n',its);

    fclose(fileID);
end

return;
%%END OF RunRSDFT FUNCTION

% function used to change the percentage complete that the progress bar shows
function updateProgressBar(listboxArray,handles,percentage,barHandle)
set(handles.OutputWindow,'String',listboxArray);
scrollPosition=max(1,numel(listboxArray)-20);
set(handles.OutputWindow,'ListboxTop',scrollPosition);

updatePercentageComplete(percentage,barHandle,handles);
