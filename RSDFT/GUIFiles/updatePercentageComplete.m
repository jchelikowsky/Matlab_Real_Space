function [] = updatePercentageComplete(percentageComplete,barHandle,handles)

set(barHandle,'YData',percentageComplete);
%axis(handles.ProgressBar,[0 1 0 .1]);
%axis(handles.ProgressBar, 'off');

drawnow();

return;