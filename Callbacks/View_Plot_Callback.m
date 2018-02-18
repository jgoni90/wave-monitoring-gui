function View_Plot_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
if value
    set(data.plotid,'Visible','on')
    if ~isempty(data.Pointer)
        set(data.Pointer,'Visible','on')
    end
    set(handles.View_PlotTransparency,'Enable','on')
else
    set(data.plotid,'Visible','off')
    if ~isempty(data.Pointer)
        set(data.Pointer,'Visible','off')
    end
    set(handles.View_PlotTransparency,'Enable','off')
end