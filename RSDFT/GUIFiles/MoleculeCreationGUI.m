function varargout = MoleculeCreationGUI(varargin)
% MOLECULECREATIONGUI M-file for MoleculeCreationGUI.fig
%      MOLECULECREATIONGUI, by itself, creates a new MOLECULECREATIONGUI or raises the existing
%      singleton*.
%
%      H = MOLECULECREATIONGUI returns the handle to a new MOLECULECREATIONGUI or the handle to
%      the existing singleton*.
%
%      MOLECULECREATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOLECULECREATIONGUI.M with the given input arguments.
%
%      MOLECULECREATIONGUI('Property','Value',...) creates a new MOLECULECREATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MoleculeCreationGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MoleculeCreationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MoleculeCreationGUI

% Last Modified by GUIDE v2.5 28-Jul-2008 09:22:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MoleculeCreationGUI_OpeningFcn, ...
    'gui_OutputFcn',  @MoleculeCreationGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MoleculeCreationGUI is made visible.
function MoleculeCreationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global CG_prec poldeg diagmeth adaptiveScheme enableChargeDensityVisualization
global OPTIMIZATIONLEVEL
global enableOSCheckForVisualization

% Choose default command line output for MoleculeCreationGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% creates the AtomsInMolecules variable and attachs it to handles
% AtomsInMolecules: an array of structs with the elements 'AtomElement' and 'Coordinate'
% unlike Atoms, each atom has its own struct
handles.AtomsInMolecule=struct('AtomElement',{},'Coordinate',{});
handles.calculationsInProgress=0;


global optWindowOpen helpWindowOpen settingsWindowOpen
optWindowOpen=0;
helpWindowOpen=0;
settingsWindowOpen=0;

guidata(hObject, handles);% needed to save AtomsInMolecule changes

% allow visualization to be rotated
rotate3d(handles.VisualizationOfMolecule);

% set GUI objects to state that reflect the values that are in the varibles
% set in the script include.m
% if statements check that values set in include.m are valid
if (CG_prec~=0 && CG_prec~=1)
    CG_prec=0;
end
if (poldeg<1)
    poldeg=1;
end
if (adaptiveScheme~=0 && adaptiveScheme~=1)
    adaptiveScheme=0;
end
if (enableChargeDensityVisualization~=0 && enableChargeDensityVisualization~=1)
    enableChargeDensityVisualization=0;
end
if (diagmeth<0 || diagmeth>3)
    diagmeth=3;
end
if (OPTIMIZATIONLEVEL~=0 && OPTIMIZATIONLEVEL~=1)
    OPTIMIZATIONLEVEL=0;
end

% disable charge density visualization if unix machine
if (enableOSCheckForVisualization~=0 && isunix()==1)
    enableChargeDensityVisualization=0;
    set(handles.visualizationCheckBox,'Enable','off');
end

% --- Outputs from this function are returned to the command line.
function varargout = MoleculeCreationGUI_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ListOfCurrentAtoms.
% function does nothing at this time
function ListOfCurrentAtoms_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
% function does nothing interesting at this time
function ListOfCurrentAtoms_CreateFcn(hObject, eventdata, handles)
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function is called when text-box is changed, selected, etc.
function XCoordinate_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of XCoordinate as text
%        str2double(get(hObject,'String')) returns contents of XCoordinate
%        as a double
xCoordinate=str2double(get(hObject,'string'));
if isnan(xCoordinate)
    errordlg('You must enter a numeric value','Bad Input','modal')
    set(hObject,'string',0);% sets value in text-box to 0
    return;
end


% --- Executes during object creation, after setting all properties.
function XCoordinate_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% initialize value of Text Box
set(hObject,'string',0);

% function is called when text-box is changed, selected, etc.
function YCoordinate_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of YCoordinate as text
%        str2double(get(hObject,'String')) returns contents of YCoordinate as a double
yCoordinate=str2double(get(hObject,'string'));
if isnan(yCoordinate)
    errordlg('You must enter a numeric value','Bad Input','modal')
    set(hObject,'string',0);% sets value in text-box to 0
    return;
end

% --- Executes during object creation, after setting all properties.
function YCoordinate_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% initialize value of Text Box
set(hObject,'string',0);


% function is called when text-box is changed, selected, etc.
function ZCoordinate_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ZCoordinate as text
%        str2double(get(hObject,'String')) returns contents of ZCoordinate as a double
zCoordinate=str2double(get(hObject,'string'));
if isnan(zCoordinate)
    errordlg('You must enter a numeric value','Bad Input','modal')
    set(hObject,'string',0);% sets value in text-box to 0
    return;
end

% --- Executes during object creation, after setting all properties.
function ZCoordinate_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% initialize value of Text Box
set(hObject,'string',0);

% --- Executes on button press in AddAtomToListButton.
function AddAtomToListButton_Callback(hObject, eventdata, handles)
if (handles.calculationsInProgress==0)
    % retrieves the contents of the coordinate text-boxes and converts strings
    % to doubles and places in an array
    coordinate=[str2double(get(handles.XCoordinate,'string')) str2double(get(handles.YCoordinate,'string')) str2double(get(handles.ZCoordinate,'string'))];
    elementOfNextAtom = handles.elementOfNextAtom;% retrieves the string representing the next atom's element

    coordinate(isnan(coordinate))=0;% replace any coordinates that are not a number with 0

    % add next atom to molecule, AtomsInMolecule is a array of structs
    AtomsInMolecule=handles.AtomsInMolecule;
    AtomsInMolecule(1,numel(AtomsInMolecule)+1)=struct('AtomElement',elementOfNextAtom,'Coordinate',coordinate);
    % save data to guidata
    handles.AtomsInMolecule=AtomsInMolecule;
    guidata(hObject, handles);

    UpdateGUIObjects(handles,AtomsInMolecule);
end


% --- Executes on button press in DeleteAtomsButton.
function DeleteAtomsButton_Callback(hObject, eventdata, handles)
if (handles.calculationsInProgress==0)
    % reset text in charge changer text box to 0
    set(handles.changeInElectronsEdit,'String','0');
    indexsOfSelectedAtoms=get(handles.ListOfCurrentAtoms,'Value');
    
    AtomsInMolecule=handles.AtomsInMolecule;

    if (numel(AtomsInMolecule)>0)
        AtomsInMolecule(indexsOfSelectedAtoms)=[];
        % save data to guidata
        handles.AtomsInMolecule=AtomsInMolecule;
        guidata(hObject, handles);

        UpdateGUIObjects(handles,AtomsInMolecule);
    end
end

% --- Executes on button press in LoadFromFileButton.
function LoadFromFileButton_Callback(hObject, eventdata, handles)
if (handles.calculationsInProgress==0)
    % reset text in charge changer text box to 0
    set(handles.changeInElectronsEdit,'String','0');
    
    % recieve the name of the file to load from uigetfile
    % will load both binary and text-based files
    defaultFileName=horzcat('.',filesep,'SavedMolecules');

    [fileName path filterIndex]=uigetfile({'*.mat' 'Binary File';'*.dat' 'Text File'},'Load Atom Configuration',defaultFileName);
    filePathName=horzcat(path,fileName);

    if (filterIndex==1)
        set(handles.changeInElectronsEdit,'String','0');
        AtomsInMolecule=struct('AtomElement',{},'Coordinate',{});

        variableNames=who('-file',filePathName);

        if (any(strcmp(variableNames,'AtomsInMolecule')))
            load(filePathName,'AtomsInMolecule');
        else
            errordlg('File does not contain correct variable','Bad File','modal')
        end

        handles.AtomsInMolecule=AtomsInMolecule;
        guidata(hObject, handles);

        UpdateGUIObjects(handles,AtomsInMolecule);
    elseif (filterIndex==2)
        set(handles.changeInElectronsEdit,'String','0');
        Atoms=struct('typ',{},'coord',{});
        [fid,message] = fopen(filePathName,'r');

        if (fid==-1)
            error(horzcat('Could not open file. ',message));
        end

        at = fscanf(fid,'%g');
        for i=1:at
            typ = fscanf(fid,'%s',1);
            readxyz = fscanf(fid,'%g');
            xyz = reshape(readxyz,3,[])';
            str=struct('typ',typ,'coord',xyz);
            Atoms(i) = str;
        end
        status = fclose(fid);
        if (status==-1)
            error('Error on closing of file');
        end
        % convert from Atoms to AtomsInMolecule
        AtomsInMolecule=struct('AtomElement',{},'Coordinate',{});

        for indexOfAtoms=1:numel(Atoms)
            for indexOfCoord=1:numel(Atoms(indexOfAtoms).coord)/3
                AtomsInMolecule(1,numel(AtomsInMolecule)+1)=struct('AtomElement',Atoms(indexOfAtoms).typ,'Coordinate',Atoms(indexOfAtoms).coord(indexOfCoord,:));
            end
        end

        handles.AtomsInMolecule=AtomsInMolecule;
        guidata(hObject, handles);

        UpdateGUIObjects(handles,AtomsInMolecule);
    end
end

% --- Executes on button press in SaveToFileButton.
function SaveToFileButton_Callback(hObject, eventdata, handles)
if (handles.calculationsInProgress==0)
    AtomsInMolecule=handles.AtomsInMolecule;
    Atoms=AtomsInMoleculeToAtomsConverter(AtomsInMolecule);
    if (numel(Atoms)>0)

        % recieve the name of the file to save to from uiputfile
        % will save both binary and text-based files
        defaultFileName=horzcat('.',filesep,'SavedMolecules');

        [fileName path filterIndex]=uiputfile({'*.mat' 'Binary File';'*.dat' 'Text File'},'Save Atom Configuration',defaultFileName);
        filePathName=horzcat(path,fileName);

        if (filterIndex==1)
            save(horzcat(filePathName,'.mat'),'AtomsInMolecule');
        elseif (filterIndex==2)
            [fid,message] = fopen(horzcat(filePathName,'.dat'),'w');

            if (fid==-1)
                error(horzcat('Could not open file. ',message));
            end

            fprintf(fid,'%g\n',numel(Atoms));
            for i=1:numel(Atoms)
                fprintf(fid,'%s\n',Atoms(i).typ);
                for indexOfCoord=1:numel(Atoms(i).coord)/3
                    fprintf(fid,'%g %g %g\n',Atoms(i).coord(indexOfCoord,1),Atoms(i).coord(indexOfCoord,2),Atoms(i).coord(indexOfCoord,3));
                end
            end
            status = fclose(fid);

            if (status==-1)
                error('Error on closing of file');
            end
        end
    else
        errordlg('You must enter at least one atom','Bad Input','modal')
    end
end

% --- Executes on button press in StartCalculationsButton.
function StartCalculationsButton_Callback(hObject, eventdata, handles)
global OPTIMIZATIONLEVEL enableChargeDensityVisualization
global enableOSCheckForVisualization

if (handles.calculationsInProgress==0)
    mexFileSuccess=1;
    if (OPTIMIZATIONLEVEL~=0)
        mexFileSuccess=TestMexFiles();
    end

    if (mexFileSuccess==1)
        handles.calculationsInProgress=1;
        guidata(hObject, handles);
        % convert to data structure format main.m uses,
        % all atoms of the same element are grouped
        AtomsInMolecule=handles.AtomsInMolecule;

        Atoms=AtomsInMoleculeToAtomsConverter(AtomsInMolecule);

        
        Z_charge=str2double(get(handles.changeInElectronsEdit,'String'));
        if (isnan(Z_charge))
            Z_charge=0;
        end
        Z_charge=-1*round(Z_charge);% * -1 to change from change in number of electrons to change in charge

        if (numel(Atoms)<=0)
            errordlg('You must enter at least one atom','Bad Input','modal')
        elseif (getNumberOfValenceElectrons(Atoms)-Z_charge<0)
            errordlg('Current System setup has less than 0 electrons','Bad Input','modal')
        else
            progressGUIhandles=ProgressGUI();
            handles.ProgressGUIhandles=progressGUIhandles;

            try
                chargeDensity=RunRSDFT(Atoms,progressGUIhandles,Z_charge);
                if (enableChargeDensityVisualization==1 && isempty(chargeDensity)==0 && (enableOSCheckForVisualization==0 || isunix()==0))
                    ChargeVisualization(chargeDensity);
                end
            catch
                exceptionID=lasterror();

                logError(exceptionID);

                warndlg({'An error occurred during execution.',horzcat('ERROR MESSAGE: ',exceptionID.message),' ','Please check your settings and if this problem persists, please contact the proper people.'},'Exception Thrown');
            end
            guidata(hObject, handles);
        end

        handles.calculationsInProgress=0;
        guidata(hObject, handles);
    else
        errordlg('Mex file(s) are missing or are incompatible with the current MATLAB version. Recompile or decrease optimization level to resolve error.','Mex File Error','modal')
    end
end

% --- Executes on button press in ExitProgramButton.
function ExitProgramButton_Callback(hObject, eventdata, handles)
close(handles.figure1);


% --- Executes on selection change in ElementSelector.
function ElementSelector_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns ElementSelector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ElementSelector
val = get(hObject,'Value');
string_list = get(hObject,'String');
handles.elementOfNextAtom= string_list{val};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ElementSelector_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% add the names of the elements to the drop-down menu
temp=getElementData(2);
set(hObject,'String',temp);
handles.elementOfNextAtom= temp{1};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function VisualizationOfMolecule_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate VisualizationOfMolecule



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global ISPROGRESSWINDOWOPEN;
% this if statement is here in the event the ProgressGUI window handle is
% not passed correctly to MoleculeCreationGUI.
if (isfield(handles,'ProgressGUIhandles'))
    otherHandle=handles.ProgressGUIhandles;

    % check to see if the window is still open
    if (ISPROGRESSWINDOWOPEN==1)
        close(otherHandle.figure1);
    end
end

global optWindowOpen helpWindowOpen settingsWindowOpen
if (optWindowOpen==1)
    tempHandle=handles.optHandles;
    close(tempHandle.figure1);
end
if (settingsWindowOpen==1)
    tempHandle=handles.settingsHandles;
    close(tempHandle.figure1);
end
if (helpWindowOpen==1)
    tempHandle=handles.helpHandles;
    close(tempHandle.figure1);
end

display 'Exiting RSDFT'

delete(hObject);


function poldegEdit_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function poldegEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poldegEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in adaptiveCheckBox.
function adaptiveCheckBox_Callback(hObject, eventdata, handles)


% --- Executes on button press in visualizationCheckBox.
function visualizationCheckBox_Callback(hObject, eventdata, handles)
global enableChargeDensityVisualization
enableChargeDensityVisualization=get(handles.visualizationCheckBox,'Value');


% --------------------------------------------------------------------
function helpMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function fileOps_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function saveMenu_Callback(hObject, eventdata, handles)
SaveToFileButton_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function loadMenu_Callback(hObject, eventdata, handles)
LoadFromFileButton_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function exitMenu_Callback(hObject, eventdata, handles)
ExitProgramButton_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function advancedMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function advancedSettingMenu_Callback(hObject, eventdata, handles)
global settingsWindowOpen
settingsWindowOpen=1;
tempHandle=AdvancedSettings;

handles.settingsHandles=tempHandle;
guidata(hObject, handles);


% --------------------------------------------------------------------
function openHelpDialog_Callback(hObject, eventdata, handles)
global helpWindowOpen
helpWindowOpen=1;
tempHandle=HelpGUI;

handles.helpHandles=tempHandle;
guidata(hObject, handles);

% --------------------------------------------------------------------
function optimizationSettingsMenu_Callback(hObject, eventdata, handles)
global optWindowOpen
optWindowOpen=1;
tempHandle=OptimizationSettings;

handles.optHandles=tempHandle;
guidata(hObject, handles);




function changeInElectronsEdit_Callback(hObject, eventdata, handles)

currentString=get(handles.changeInElectronsEdit,'String');
currentValue=str2double(currentString);
currentValue=round(currentValue);

roundedString=sprintf('%d',currentValue);
set(handles.changeInElectronsEdit,'String',roundedString);

% test for valid input
if (isnan(currentValue)==1)
    set(handles.changeInElectronsEdit,'String','0');
end


% --- Executes during object creation, after setting all properties.
function changeInElectronsEdit_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addElectron.
function addElectron_Callback(hObject, eventdata, handles)

% test if another electron can be added
currentString=get(handles.changeInElectronsEdit,'String');
currentValue=str2double(currentString);
if (isnan(currentValue)==0)
    newValue=currentValue+1;
    newString=sprintf('%d',newValue);
    set(handles.changeInElectronsEdit,'String',newString);
end




% --- Executes on button press in removeElectron.
function removeElectron_Callback(hObject, eventdata, handles)
AtomsInMolecule=handles.AtomsInMolecule;
Atoms=AtomsInMoleculeToAtomsConverter(AtomsInMolecule);

currentString=get(handles.changeInElectronsEdit,'String');
currentValue=str2double(currentString);
% test if an electron can be removed
if (isnan(currentValue)==0 && getNumberOfValenceElectrons(Atoms)+currentValue>0)
        newValue=currentValue-1;
        newString=sprintf('%d',newValue);
        set(handles.changeInElectronsEdit,'String',newString);
end

