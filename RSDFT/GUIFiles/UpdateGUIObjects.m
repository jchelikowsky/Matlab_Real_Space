function success=UpdateGUIObjects(handles,AtomsInMolecule)
% called in MoleculeCreationGUI when ever the GUI Objects need to change
% what they are displaying.  Updates ListOfCurrentAtoms and the plot of the
% Atoms' coordinates

success=1;

% update the strings in the List to reflect changes in AtomsInMolecule
cellArrayOfStrings=convertStructToString(AtomsInMolecule);
set(handles.ListOfCurrentAtoms,'String',cellArrayOfStrings,'Value',1);       
set(handles.ListOfCurrentAtoms,'Max',numel(cellArrayOfStrings));
% Max-Min is the number of items that can be concurrently selected.
% Min defaults to 1


% graph each atom to the plot with the element name next to the data
% point
elementGraphicsData=getElementData('graphics');
elementNameData=getElementData('names');
    
cla;
hold on;
for box=1:numel(AtomsInMolecule);
    coordinates=AtomsInMolecule(box).Coordinate;
    plot3(coordinates(1),coordinates(2),coordinates(3),elementGraphicsData{strcmp(elementNameData,AtomsInMolecule(box).AtomElement)});
    text(coordinates(1),coordinates(2),coordinates(3),elementNameData{strcmp(elementNameData,AtomsInMolecule(box).AtomElement)});
end    
hold off;



return