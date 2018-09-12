function varargout = AdvancedSettings(varargin)
% ADVANCEDSETTINGS M-file for AdvancedSettings.fig
%      ADVANCEDSETTINGS, by itself, creates a new ADVANCEDSETTINGS or raises the existing
%      singleton*.
%
%      H = ADVANCEDSETTINGS returns the handle to a new ADVANCEDSETTINGS or the handle to
%      the existing singleton*.
%
%      ADVANCEDSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDSETTINGS.M with the given input arguments.
%
%      ADVANCEDSETTINGS('Property','Value',...) creates a new ADVANCEDSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdvancedSettings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdvancedSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdvancedSettings

% Last Modified by GUIDE v2.5 13-Aug-2009 11:35:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdvancedSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @AdvancedSettings_OutputFcn, ...
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


% --- Executes just before AdvancedSettings is made visible.
function AdvancedSettings_OpeningFcn(hObject, eventdata, handles, varargin)
global CG_prec poldeg diagmeth adaptiveScheme
global fd_order tol maxits Fermi_temp

% Choose default command line output for AdvancedSettings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set GUI objects to state that reflect the values that are in the varibles
% set in the script settings.m
% if statements check that values set in settings.m are valid
if (CG_prec~=0 && CG_prec~=1)
    CG_prec=0;
end
if (poldeg<1)
    poldeg=1;
end
if (adaptiveScheme~=0 && adaptiveScheme~=1)
    adaptiveScheme=0;
end
if (diagmeth<0 || diagmeth>3)
    diagmeth=3;
end
fd_order=round(fd_order);
if (fd_order<1)
    fd_order=1;
end
if (tol<0)
    tol=0.000001; 
end
maxits=round(maxits);
if (maxits<1)
    maxits=1;
end    
if (Fermi_temp<1)
   Fermi_temp=1; 
end
    
set(handles.PreconditionCG,'Value',CG_prec);
polString=sprintf('%d',poldeg);
set(handles.degreeTextBox,'String',polString);
set(handles.schemeCheckBox,'Value',adaptiveScheme);
set(handles.diagMethodhSelect,'Value',diagmeth+1);

fdString=sprintf('%d',fd_order);
set(handles.fdorderEdit,'String',fdString);
tolString=sprintf('%f',tol);
set(handles.tolerenceEdit,'String',tolString);
maxitsString=sprintf('%d',maxits);
set(handles.maxitsEdit,'String',maxitsString);
fermiString=sprintf('%f',Fermi_temp);
set(handles.fermiTempEdit,'String',fermiString);

% --- Outputs from this function are returned to the command line.
function varargout = AdvancedSettings_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles;


% --- Executes on button press in PreconditionCG.
function PreconditionCG_Callback(hObject, eventdata, handles)
global CG_prec
CG_prec=get(handles.PreconditionCG,'Value');


% --- Executes on button press in schemeCheckBox.
function schemeCheckBox_Callback(hObject, eventdata, handles)
global adaptiveScheme
adaptiveScheme=get(handles.schemeCheckBox,'Value');


% --- Executes on selection change in diagMethodhSelect.
function diagMethodhSelect_Callback(hObject, eventdata, handles)
global diagmeth
diagmeth=get(handles.diagMethodhSelect,'Value')-1;


% --- Executes during object creation, after setting all properties.
function diagMethodhSelect_CreateFcn(hObject, eventdata, handles)


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function degreeTextBox_Callback(hObject, eventdata, handles)
global poldeg
polynomialDegree=str2double(get(handles.degreeTextBox,'String'));
if (polynomialDegree<1)
    polynomialDegree=1;
end
poldeg=polynomialDegree;


% --- Executes during object creation, after setting all properties.
function degreeTextBox_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in defaultButton.
function defaultButton_Callback(hObject, eventdata, handles)
global OPTIMIZATIONLEVEL enableChargeDensityVisualization

optimizeCopy=OPTIMIZATIONLEVEL;
visualizationCopy=enableChargeDensityVisualization;


include;
% set GUI objects to state that reflect the values that are in the varibles
% set in the script settings.m
% if statements check that values set in settings.m are valid
if (CG_prec~=0 && CG_prec~=1)
    CG_prec=0;
end
if (poldeg<1)
    poldeg=1;
end
if (adaptiveScheme~=0 && adaptiveScheme~=1)
    adaptiveScheme=0;
end
if (diagmeth<0 || diagmeth>3)
    diagmeth=3;
end
fd_order=round(fd_order);
if (fd_order<1)
    fd_order=1;
end
if (tol<0)
    tol=0.000001; 
end
maxits=round(maxits);
if (maxits<1)
    maxits=1;
end    
if (Fermi_temp<1)
   Fermi_temp=1; 
end
set(handles.PreconditionCG,'Value',CG_prec);
polString=sprintf('%d',poldeg);
set(handles.degreeTextBox,'String',polString);
set(handles.schemeCheckBox,'Value',adaptiveScheme);
set(handles.diagMethodhSelect,'Value',diagmeth+1);

fdString=sprintf('%d',fd_order);
set(handles.fdorderEdit,'String',fdString);
tolString=sprintf('%f',tol);
set(handles.tolerenceEdit,'String',tolString);
maxitsString=sprintf('%d',maxits);
set(handles.maxitsEdit,'String',maxitsString);
fermiString=sprintf('%f',Fermi_temp);
set(handles.fermiTempEdit,'String',fermiString);

OPTIMIZATIONLEVEL=optimizeCopy;
enableChargeDensityVisualization=visualizationCopy;




function fdorderEdit_Callback(hObject, eventdata, handles)
global fd_order
newValue=str2double(get(handles.fdorderEdit,'String'));
newValue=round(newValue);
if (newValue<1)
    newValue=1;
end
fd_order=newValue;
    

% --- Executes during object creation, after setting all properties.
function fdorderEdit_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tolerenceEdit_Callback(hObject, eventdata, handles)
global tol
newValue=str2double(get(handles.tolerenceEdit,'String'));
if (newValue<0)
    newValue=.000001;
end
tol=newValue;


% --- Executes during object creation, after setting all properties.
function tolerenceEdit_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxitsEdit_Callback(hObject, eventdata, handles)
global maxits
newValue=str2double(get(handles.maxitsEdit,'String'));
newValue=round(newValue);
if (newValue<1)
    newValue=1;
end
maxits=newValue;

% --- Executes during object creation, after setting all properties.
function maxitsEdit_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fermiTempEdit_Callback(hObject, eventdata, handles)
global Fermi_temp
newValue=str2double(get(handles.fermiTempEdit,'String'));
if (newValue<1)
    newValue=1;
end
Fermi_temp=newValue;


% --- Executes during object creation, after setting all properties.
function fermiTempEdit_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global settingsWindowOpen
settingsWindowOpen=0;

% Hint: delete(hObject) closes the figure
delete(hObject);
