function outputArray = getElementData(dataDesired)
% This function returns data about elements for use in GUIs
% the only input is used to specify which data is desired, can be a string
% or a number
% if an element is added to the program, this is one of the files that must
% be modified


elementNames={'Al';'Ar';'B';'Be';'C';'Cl';'F';'H';'He';'Li';'Mg';'N';'Na';'Ne';'O';'P';'S';'Si'};
% elementNames above, elementGraphics below
elementGraphics={'r*';'g*';'b*';'r+';'g+';'b+';'rs';'gs';'bs';'r^';'g^';'b^';'ro';'go';'bo';'rp';'gp';'bp'};

% add data and elseif's as needed
elementData=struct('Names',elementNames,'Graphics',elementGraphics);

if ((ischar(dataDesired) && strcmpi(dataDesired,'all')) || (isa(dataDesired,'numeric') && dataDesired==1))
    outputArray=elementData;
elseif ((ischar(dataDesired) && strcmpi(dataDesired,'names')) || (isa(dataDesired,'numeric') && dataDesired==2))
    outputArray=elementNames;
elseif ((ischar(dataDesired) && strcmpi(dataDesired,'graphics')) || (isa(dataDesired,'numeric') && dataDesired==3))
    outputArray=elementGraphics;
else
    outputArray=[];
end

return;