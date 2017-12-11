function varargout = HelpGUI(varargin)
% HELPGUI M-file for HelpGUI.fig
%      HELPGUI, by itself, creates a new HELPGUI or raises the existing
%      singleton*.
%
%      H = HELPGUI returns the handle to a new HELPGUI or the handle to
%      the existing singleton*.
%
%      HELPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELPGUI.M with the given input arguments.
%
%      HELPGUI('Property','Value',...) creates a new HELPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HelpGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HelpGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HelpGUI

% Last Modified by GUIDE v2.5 26-Jul-2008 13:13:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HelpGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HelpGUI_OutputFcn, ...
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


% --- Executes just before HelpGUI is made visible.
function HelpGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% set helpInfo text list to show first topic on open
topicString=get(handles.topicList,'String');
helpText=loadHelpText(topicString{1});
set(handles.helpInfo,'String',helpText);


% --- Outputs from this function are returned to the command line.
function varargout = HelpGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles;


% --- Executes on selection change in topicList.
function topicList_Callback(hObject, eventdata, handles)
topicString=get(handles.topicList,'String');
index=get(handles.topicList,'Value');
helpText=loadHelpText(topicString{index});
set(handles.helpInfo,'String',helpText);


% --- Executes during object creation, after setting all properties.
function topicList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in helpInfo.
function helpInfo_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function helpInfo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global helpWindowOpen
helpWindowOpen=0;

% Hint: delete(hObject) closes the figure
delete(hObject);


