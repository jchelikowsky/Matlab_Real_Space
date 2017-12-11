function varargout = OptimizationSettings(varargin)
% OPTIMIZATIONSETTINGS M-file for OptimizationSettings.fig
%      OPTIMIZATIONSETTINGS, by itself, creates a new OPTIMIZATIONSETTINGS or raises the existing
%      singleton*.
%
%      H = OPTIMIZATIONSETTINGS returns the handle to a new OPTIMIZATIONSETTINGS or the handle to
%      the existing singleton*.
%
%      OPTIMIZATIONSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTIMIZATIONSETTINGS.M with the given input arguments.
%
%      OPTIMIZATIONSETTINGS('Property','Value',...) creates a new OPTIMIZATIONSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OptimizationSettings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OptimizationSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OptimizationSettings

% Last Modified by GUIDE v2.5 26-Jul-2008 13:14:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OptimizationSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @OptimizationSettings_OutputFcn, ...
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


% --- Executes just before OptimizationSettings is made visible.
function OptimizationSettings_OpeningFcn(hObject, eventdata, handles, varargin)
global OPTIMIZATIONLEVEL

% Choose default command line output for OptimizationSettings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if (OPTIMIZATIONLEVEL~=0 && OPTIMIZATIONLEVEL~=1)
   OPTIMIZATIONLEVEL=0; 
end 
set(handles.optimizationLevelSelect,'Value',OPTIMIZATIONLEVEL+1);


% --- Outputs from this function are returned to the command line.
function varargout = OptimizationSettings_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes on selection change in optimizationLevelSelect.
function optimizationLevelSelect_Callback(hObject, eventdata, handles)
global OPTIMIZATIONLEVEL;
value=get(hObject,'Value');

if (value==1)
    OPTIMIZATIONLEVEL=0;
elseif (value==2)
    OPTIMIZATIONLEVEL=1;
else
    OPTIMIZATIONLEVEL=0;
end


% --- Executes during object creation, after setting all properties.
function optimizationLevelSelect_CreateFcn(hObject, eventdata, handles)


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in compileButton.
function compileButton_Callback(hObject, eventdata, handles)
compileMexFiles(1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global optWindowOpen
optWindowOpen=0;

% Hint: delete(hObject) closes the figure
delete(hObject);


