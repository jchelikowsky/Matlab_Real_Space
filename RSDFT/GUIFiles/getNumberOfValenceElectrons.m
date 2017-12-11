function [zelec] = getNumberOfValenceElectrons(Atoms)

if (isempty(Atoms))
   zelec=0;
   return;
end    

n_atom=zeros(1,numel(Atoms));
for tempIndex=1:numel(Atoms)
    n_atom(tempIndex)=numel(Atoms(tempIndex).coord)/3;
end

elem    = importdata('elements_new.csv');
N_elements=size(elem.data);
%%%%
%%%%
% This is the number of species that need to be looked at
N_types = length(Atoms);


zelec=0.;
for at_typ=1:N_types
    % This gets the first atom and assigns the atomic symbol to typ
    typ = Atoms(at_typ).typ;
    % This iterates through elem looking for the matching element data

    %%% Loop over number of atoms
    for i = 1: N_elements
        if strcmp(typ,elem.textdata{i})
            index=i;
            Z=elem.data(index,2)*n_atom(at_typ);
            zelec=zelec+Z;
        end
    end
end

return;