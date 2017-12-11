% Spline data for atoms
% data is loaded from the binary file splineData.mat
% file should be compatable with 7.* versions of MATLAB
% This is one of the files that must be modified if an atom is
% added to the program

nAtoms = 0;
% H data
nAtoms = nAtoms + 1;
atom = 'H';
load('splineData.mat','dataH');
str = struct('atom', atom, 'data', dataH);
AtomFuncData(nAtoms) = str;

% He data
nAtoms = nAtoms + 1;
atom = 'He';
load('splineData.mat','dataHe');
str = struct('atom', atom, 'data', dataHe);
AtomFuncData(nAtoms) = str;

% Li data
nAtoms = nAtoms + 1;
atom = 'Li';
load('splineData.mat','dataLi');
str = struct('atom', atom, 'data', dataLi);
AtomFuncData(nAtoms) = str;

% Be data
nAtoms = nAtoms + 1;
atom = 'Be';
load('splineData.mat','dataBe');
str = struct('atom', atom, 'data', dataBe);
AtomFuncData(nAtoms) = str;

% B data
nAtoms = nAtoms + 1;
atom = 'B';
load('splineData.mat','dataB');
str = struct('atom', atom, 'data', dataB);
AtomFuncData(nAtoms) = str;

% C data
nAtoms = nAtoms + 1;
atom = 'C';
load('splineData.mat','dataC');
str = struct('atom', atom, 'data', dataC);
AtomFuncData(nAtoms) = str;

% N data
nAtoms = nAtoms + 1;
atom = 'N';
load('splineData.mat','dataN');
str = struct('atom', atom, 'data', dataN);
AtomFuncData(nAtoms) = str;

% O data
nAtoms = nAtoms + 1;
atom = 'O';
load('splineData.mat','dataO');
str = struct('atom', atom, 'data', dataO);
AtomFuncData(nAtoms) = str;

% F data
nAtoms = nAtoms + 1;
atom = 'F';
load('splineData.mat','dataF');
str = struct('atom', atom, 'data', dataF);
AtomFuncData(nAtoms) = str;

% Ne data
nAtoms = nAtoms + 1;
atom = 'Ne';
load('splineData.mat','dataNe');
str = struct('atom', atom, 'data', dataNe);
AtomFuncData(nAtoms) = str;

% Na data
nAtoms = nAtoms + 1;
atom = 'Na';
load('splineData.mat','dataNa');
str = struct('atom', atom, 'data', dataNa);
AtomFuncData(nAtoms) = str;

% Mg data
nAtoms = nAtoms + 1;
atom = 'Mg';
load('splineData.mat','dataMg');
str = struct('atom', atom, 'data', dataMg);
AtomFuncData(nAtoms) = str;

% Al data
nAtoms = nAtoms + 1;
atom = 'Al';
load('splineData.mat','dataAl');
str = struct('atom', atom, 'data', dataAl);
AtomFuncData(nAtoms) = str;

% Si data
nAtoms = nAtoms + 1;
atom = 'Si';
load('splineData.mat','dataSi');
str = struct('atom', atom, 'data', dataSi);
AtomFuncData(nAtoms) = str;

% P data
nAtoms = nAtoms + 1;
atom = 'P';
load('splineData.mat','dataP');
str = struct('atom', atom, 'data', dataP);
AtomFuncData(nAtoms) = str;

% S data
nAtoms = nAtoms + 1;
atom = 'S';
load('splineData.mat','dataS');
str = struct('atom', atom, 'data', dataS);
AtomFuncData(nAtoms) = str;

% Cl data
nAtoms = nAtoms + 1;
atom = 'Cl';
load('splineData.mat','dataCl');
str = struct('atom', atom, 'data', dataCl);
AtomFuncData(nAtoms) = str;

% Ar data
nAtoms = nAtoms + 1;
atom = 'Ar';
load('splineData.mat','dataAr');
str = struct('atom', atom, 'data', dataAr);
AtomFuncData(nAtoms) = str;


%add new elements at this line


% List of data name:
data_list = [ {'radius'}, {'charge'}, {'hartree'}, {'pot_P'}, {'pot_S'}, {'wfn_P'}, {'wfn_S'} ];