function varargout = ProgressGUI(varargin)
% PROGRESSGUI M-file for ProgressGUI.fig
%      PROGRESSGUI, by itself, creates a new PROGRESSGUI or raises the existing
%      singleton*.
%
%      H = PROGRESSGUI returns the handle to a new PROGRESSGUI or the handle to
%      the existing singleton*.
%
%      PROGRESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROGRESSGUI.M with the given input arguments.
%
%      PROGRESSGUI('Property','Value',...) creates a new PROGRESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProgressGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProgressGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProgressGUI

% Last Modified by GUIDE v2.5 12-Jun-2008 20:28:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProgressGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ProgressGUI_OutputFcn, ...
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


% --- Executes just before ProgressGUI is made visible.
function ProgressGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ProgressGUI

global ISPROGRESSWINDOWOPEN;
ISPROGRESSWINDOWOPEN=1;

% sets up the axis
axis([0 1 0 .1])
axis off;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProgressGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProgressGUI_OutputFcn(hObject, eventdata, handles) 

% return the handles structure to calling function
% gives calling function access to all variables in handles
varargout{1} = handles;

% --- Executes on selection change in OutputWindow.
function OutputWindow_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns OutputWindow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OutputWindow


% --- Executes during object creation, after setting all properties.
function OutputWindow_CreateFcn(hObject, eventdata, handles)
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global ISPROGRESSWINDOWOPEN;
ISPROGRESSWINDOWOPEN=0;

delete(hObject);



% --- Executes during object creation, after setting all properties.
function ProgressBar_CreateFcn(hObject, eventdata, handles)
% Hint: place code in OpeningFcn to populate ProgressBar

