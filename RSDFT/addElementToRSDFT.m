function [] = addElementToRSDFT()

disp('This function is is designed to automate the process of');
disp('adding elements in RSDFT');
disp('It is suggested that all the files that will be being changed');
disp('be backed-up');
disp('')

%% Input element specific data
elementName=input('What is the abbreviation of the new element? ','s');
atomicNumber=input('What is its atomic number? ');
valenceElectrons=input('How many valence Electrons does it have? ');
zValue=input('What is its Z Value? ');
gridSize=input('What is the recommended grid size for this element? ');
radius=input('What is the elements''s radius? ');
coreSize=input('What is the core size? ');

plotString=input('Input the string used to control the symbol used on the plot. i.e. r*: ','s');

%% Find file with spline data and copy to directory
disp('Now find the script that contains the spline data for this element');
disp('the variable in the script is assumed to be of the form, data*, where the *');
disp('is the element''s abbreviation');

currentDirectory=fileparts(which('addElementToRSDFT.m'));

[fileName,pathName]=uigetfile('*.m','Select the Script with the spline data');

absoluteFileName=horzcat(pathName,fileName);

destination=horzcat(currentDirectory,filesep,'SPLINES',filesep,'Spline Data Backup',filesep,fileName);

copyfile(absoluteFileName,destination);

variableName=horzcat('data',elementName);


%% Make changes to SplineDataCreator.m
filenameOfSplineDataCreator=horzcat(currentDirectory,filesep,'SPLINES',filesep,'Spline Data Backup',filesep,'SplineDataCreator.m');


try
    cellArrayOfMFile=textread(filenameOfSplineDataCreator,'%s','delimiter', '\n','whitespace', '');
catch
    exceptionID=lasterror();
   
    warndlg({'An error occurred while trying to read SplineDataCreator.m. Element Adding FAILED',horzcat('ERROR MESSAGE: ',exceptionID.message)},'Exception Thrown');  
    return;
end    


indexOfOnePastLastScript=find(strcmp('%call script to load data into a variable',cellArrayOfMFile));

scriptName=fileName(1:end-2);

% insert a row in the m file
cellArrayOfMFile(indexOfOnePastLastScript+1:end+1)=cellArrayOfMFile(indexOfOnePastLastScript:end);
cellArrayOfMFile(indexOfOnePastLastScript)=cellstr(horzcat(scriptName,';'));

indexOnePastSave=find(strcmp('%call save with the file splineData.mat and the name of the variable,',cellArrayOfMFile));

saveStatement=cellArrayOfMFile{indexOnePastSave-1};
desiredText=',''-append''';
indexInSaveStatement=findstr(desiredText,saveStatement);

% copy end of save statement and add new variable to statement
saveStatement(indexInSaveStatement+numel(variableName)+3:end+numel(variableName)+3)=saveStatement(indexInSaveStatement:end);

charactersInVariableStatement=numel(horzcat(',''',variableName,''''));
saveStatement(indexInSaveStatement:indexInSaveStatement+charactersInVariableStatement-1)=horzcat(',''',variableName,'''');

cellArrayOfMFile(indexOnePastSave-1)=cellstr(saveStatement);

% save changes to splineDataCreator
fileHandle=fopen(filenameOfSplineDataCreator,'w');
for i=1:numel(cellArrayOfMFile)
    fprintf(fileHandle,'%s\n',cellArrayOfMFile{i});
end
fclose(fileHandle);


%%  Make and move splineData.mat
% need to use p file because MATLAB will not allow an M file to
% be called after is has been modified
path(path,horzcat(currentDirectory,filesep,'SPLINES',filesep,'Spline Data Backup'));
pcode(filenameOfSplineDataCreator);
SplineDataCreator;% saves variables to a file in pwd
delete('SplineDataCreator.p');
rmpath(horzcat(currentDirectory,filesep,'SPLINES',filesep,'Spline Data Backup'));

source=horzcat(currentDirectory,filesep,'splineData.mat');
destination=horzcat(currentDirectory,filesep,'IonicPotentialFindingFiles',filesep,'splineData.mat');
movefile(source,destination);
%% modify splineData.m
cellArrayOfMFile={};
fileNameOfsplineData=horzcat(currentDirectory,filesep,'IonicPotentialFindingFiles',filesep,'splineData.m');

try
    cellArrayOfMFile=textread(fileNameOfsplineData,'%s','delimiter', '\n','whitespace', '');
catch
    exceptionID=lasterror();
   
    warndlg({'An error occurred while trying to read splineData.m. Element Adding FAILED',horzcat('ERROR MESSAGE: ',exceptionID.message)},'Exception Thrown');      
end    

indexOfLineInsplineData=find(strcmp('%add new elements at this line',cellArrayOfMFile));
cellArrayOfMFile(indexOfLineInsplineData+7:end+7)=cellArrayOfMFile(indexOfLineInsplineData:end);

cellArrayOfMFile(indexOfLineInsplineData)=cellstr(horzcat('% ',elementName, ' data'));
cellArrayOfMFile(indexOfLineInsplineData+1)=cellstr('nAtoms = nAtoms + 1;');
cellArrayOfMFile(indexOfLineInsplineData+2)=cellstr(horzcat('atom = ','''',elementName,''';'));
cellArrayOfMFile(indexOfLineInsplineData+3)=cellstr(horzcat('load(''splineData.mat'',','''',elementName,''');'));
cellArrayOfMFile(indexOfLineInsplineData+4)=cellstr(horzcat('str = struct(''atom'', atom, ''data'',', variableName,');'));
cellArrayOfMFile(indexOfLineInsplineData+5)=cellstr('AtomFuncData(nAtoms) = str;');


% save changes to splineData.m
fileHandle=fopen(fileNameOfsplineData,'w');
for i=1:numel(cellArrayOfMFile)
    fprintf(fileHandle,'%s\n',cellArrayOfMFile{i});
end
fclose(fileHandle);

%% modify element_new.csv
fileHandle=fopen('elements_new.csv','a');
fprintf(fileHandle,'%s,%d,%d,%f,%f,%f,%f\n',elementName,atomicNumber,valenceElectrons,zValue,gridSize,radius,coreSize);
fclose(fileHandle);

%% modify getElementData.m
cellArrayOfMFile={};
fileNameOfgetElementData=horzcat('.',filesep,'GUIFiles',filesep,'getElementData.m');

try
    cellArrayOfMFile=textread(fileNameOfgetElementData,'%s','delimiter', '\n','whitespace', '');
catch
    exceptionID=lasterror();
   
    warndlg({'An error occurred while trying to read getElementData.m. Element Adding FAILED',horzcat('ERROR MESSAGE: ',exceptionID.message)},'Exception Thrown');      
end    

indexOfElementData=find(strcmp('% elementNames above, elementGraphics below',cellArrayOfMFile));
elementNameStatement=cellArrayOfMFile{indexOfElementData-1};
elementGraphicsStatement=cellArrayOfMFile{indexOfElementData+1};

indexInName=findstr('}',elementNameStatement);
indexInGraphics=findstr('}',elementGraphicsStatement);

elementNameStatement(indexInName+3+numel(elementName):end+3+numel(elementName))=elementNameStatement(indexInName:end);
elementNameStatement(indexInName:indexInName+2+numel(elementName))=horzcat(';''',elementName,'''');
cellArrayOfMFile(indexOfElementData-1)=cellstr(elementNameStatement);

elementGraphicsStatement(indexInGraphics+3+numel(plotString):end+3+numel(plotString))=elementGraphicsStatement(indexInGraphics:end);
elementGraphicsStatement(indexInGraphics:indexInGraphics+2+numel(plotString))=horzcat(';''',plotString,'''');
cellArrayOfMFile(indexOfElementData+1)=cellstr(elementGraphicsStatement);

% save changes to splineData.m
fileHandle=fopen(fileNameOfgetElementData,'w');
for i=1:numel(cellArrayOfMFile)
    fprintf(fileHandle,'%s\n',cellArrayOfMFile{i});
end
fclose(fileHandle);

