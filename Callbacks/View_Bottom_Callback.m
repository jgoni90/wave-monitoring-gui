function View_Bottom_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
if value
    set(data.bottomid,'Visible','on')
    set(handles.View_BotTransparency,'Enable','on')
else
    set(data.bottomid,'Visible','off')
    set(handles.View_BotTransparency,'Enable','off')
end

guidata(hObject,data);