function View_Axes_Callback(hObject,eventdata,handles)

value = get(hObject,'Value');
if value
    set(handles.axesPlot1,'Visible','on')
    set(handles.axesPlot2,'Visible','on')
    set(handles.axesPlot3,'Visible','on')
    set(handles.axesPlot4,'Visible','on')
else
    set(handles.axesPlot1,'Visible','off')
    set(handles.axesPlot2,'Visible','off')
    set(handles.axesPlot3,'Visible','off')
    set(handles.axesPlot4,'Visible','off')
end