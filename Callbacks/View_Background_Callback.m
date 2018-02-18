function View_Background_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');

if value
    set(data.backgroundid,'Visible','on')
    set(handles.View_BGTransparency,'Enable','on')
else
    set(data.backgroundid,'Visible','off')
    set(handles.View_BGTransparency,'Enable','off')
end