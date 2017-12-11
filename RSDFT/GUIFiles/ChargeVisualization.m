function varargout = ChargeVisualization(varargin)
% CHARGEVISUALIZATION M-file for ChargeVisualization.fig
%      CHARGEVISUALIZATION, by itself, creates a new CHARGEVISUALIZATION or raises the existing
%      singleton*.
%
%      H = CHARGEVISUALIZATION returns the handle to a new CHARGEVISUALIZATION or the handle to
%      the existing singleton*.
%
%      CHARGEVISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHARGEVISUALIZATION.M with the given input arguments.
%
%      CHARGEVISUALIZATION('Property','Value',...) creates a new CHARGEVISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChargeVisualization_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChargeVisualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChargeVisualization

% Last Modified by GUIDE v2.5 02-Jul-2008 16:39:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChargeVisualization_OpeningFcn, ...
                   'gui_OutputFcn',  @ChargeVisualization_OutputFcn, ...
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


% --- Executes just before ChargeVisualization is made visible.
function ChargeVisualization_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ChargeVisualization
handles.output = hObject;

handles.chargeDensity=varargin{1};

% Update handles structure
guidata(hObject, handles);

set(handles.elevationSlider,'value',.5);
updateSurf(.5,handles.chargeDensity,handles);
rotate3d(handles.chargeVisualization,'on');
view(handles.chargeVisualization,3);
% UIWAIT makes ChargeVisualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.figure1,'Position',[5,50,140,44]);



% --- Outputs from this function are returned to the command line.
function varargout = ChargeVisualization_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on slider movement.
function elevationSlider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
position=get(handles.elevationSlider,'value');
updateSurf(position,handles.chargeDensity,handles);


% --- Executes during object creation, after setting all properties.
function elevationSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function updateSurf(position,data,handle)
global enableOSCheckForVisualization
if (enableOSCheckForVisualization==0 || isunix()==0)
    % convert from the slider position to an index in the data
    yIndex=round(position*size(data,2));
    if (yIndex<1)
        yIndex=1;
    elseif (yIndex>size(data,2))
        yIndex=size(data,2);
    end


    % get slice and save view window details so when
    % redrawing, the view does not change
    chargeDensitySlice=squeeze(data(:,yIndex,:));
    [az, el]=view;

    maxValue=max(max(max(data)));
    minValue=min(min(min(data)));

    zoomValue=get(handle.zoomSlider,'Value');

    diff=maxValue-minValue;

    if (zoomValue==1)
        zoomValue=0.999;
    end

    % interp2 smooths the data by adding points
    yInterpolation=interp2(chargeDensitySlice,1);

    surf(handle.chargeVisualization,yInterpolation);
    axis(handle.chargeVisualization,[0, size(yInterpolation,1), 0, size(yInterpolation,2), minValue, maxValue-zoomValue*diff, minValue, maxValue-zoomValue*diff]);
    view(handle.chargeVisualization,az,el);
    shading(handle.chargeVisualization,'interp');
else
    warndlg({'Visualization is disabled for unix machines due to it sometimes crashing the window manager'},'Visualization Disabled');
end

% --- Executes during object creation, after setting all properties.
function chargeVisualization_CreateFcn(hObject, eventdata, handles)




% --- Executes on slider movement.
function zoomSlider_Callback(hObject, eventdata, handles)
position=get(handles.elevationSlider,'value');
updateSurf(position,handles.chargeDensity,handles);

% --- Executes during object creation, after setting all properties.
function zoomSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


