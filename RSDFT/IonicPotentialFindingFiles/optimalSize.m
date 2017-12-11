function [imax]=optimalSize(Atoms, h)

imax = 0;

% This is just importing the information about each element
elem    = importdata('elements_new.csv');

% This is the number of atoms that need to be looked at
N_types = length(Atoms);

%% Loop over atom types
for at_typ=1:N_types
  % This gets the first atom and assigns the atomic symbol to typ
  typ = Atoms(at_typ).typ;
  % This iterates through elem looking for the matching element data
  for i = 1:length(elem.data)
    if strcmp(typ,elem.textdata{i})
      index=i;
      break
    end
  end
  % xyz corresponds to the atom's position within the Domain
  xyz=Atoms(at_typ).coord;
  % natoms is just the number of typ atoms that are in the system
  natoms = size(xyz,1);
  % Rzero is the core size of the atom
  Rzero = elem.data(index,6);
  span  = round(Rzero/h);
  imax  = imax + (4*span^3)^2 * natoms;
end % end of at_typ for loop
imax = round(imax);
end % end of function
