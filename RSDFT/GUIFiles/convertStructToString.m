function stringOfAtomData=convertStructToString(atomsStruct)
% takes in a a structure with two fields, Coordinate and AtomElement
% converts the structure to a cell array 
% function is used to generate the strings used in the list of atoms GUI
% object in MoleculeCreationGUI
lengthOfStructure=numel(atomsStruct);

stringOfAtomData=cell(lengthOfStructure,1);

for i=1:lengthOfStructure
    coordinate=atomsStruct(i).Coordinate;
    stringOfAtomData(i)=cellstr(strcat(atomsStruct(i).AtomElement,': ',...
            num2str(coordinate(1)),', ',num2str(coordinate(2)),', ',num2str(coordinate(3))));
end
return