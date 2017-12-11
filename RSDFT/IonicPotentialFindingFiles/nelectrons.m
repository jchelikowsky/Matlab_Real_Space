function [Nelec] = nelectrons(Atoms)
%% obtains the total 
%% number of valence electrons 
%% 
Nelec = 0; 
N_types = length(Atoms);
for at_typ=1:N_types
    atom = Atoms(at_typ);
    typ = atom.typ;
    elem=importdata('elements_new.csv');
    for i=1:length(elem.data)
        if strcmp(typ,elem.textdata{i})
            index=i;
            break
        end
    end
    val = elem.data(index,2) ;
    Nelec = Nelec + size(atom.coord,1)*val;
end
    
