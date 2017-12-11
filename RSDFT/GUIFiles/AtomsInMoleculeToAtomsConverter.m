function Atoms=AtomsInMoleculeToAtomsConverter(AtomsInMolecule)
% Takes a structure that holds each atom idividually and
% converts it to a structure where atoms of the same element are grouped


% listOfElements should be a cell array of strings where the strings
% are the names of the elements, ie H, He, Li, etc.
listOfElements=getElementData('names');

% add each element to the Atoms structure
Atoms=struct('typ',{},'coord',{});
for index=numel(listOfElements):-1:1
   Atoms(index).typ= listOfElements{index};
end

% take each atom in AtomsInMolecule and place it in the appropriate
% index in Atoms
for index=1:numel(AtomsInMolecule)
    for indexOfAtoms=1:numel(Atoms)
       if (strcmp(AtomsInMolecule(index).AtomElement,Atoms(indexOfAtoms).typ))
             Atoms(indexOfAtoms).coord=[Atoms(indexOfAtoms).coord;AtomsInMolecule(index).Coordinate];
             break;% exit the inner loop because only one should match
       end    
    end    
end

%remove the elements in Atoms that have no coordinates
indexsOfElementsToRemove=[];
for index=1:numel(Atoms)
    if (isempty(Atoms(index).coord))
        indexsOfElementsToRemove(numel(indexsOfElementsToRemove)+1)=index;
    end    
end
Atoms(indexsOfElementsToRemove)=[];

return;