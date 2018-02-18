function View_Mesh_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
if value
%     set(data.plotid,'EdgeAlpha',1,'FaceColor','w')
    set(data.meshid,'Visible','on')
    set(handles.View_MeshTransparency,'Enable','on')
else
%     set(data.plotid,'EdgeAlpha',0,'FaceColor','interp')
    set(data.meshid,'Visible','off')
    set(handles.View_MeshTransparency,'Enable','off')
end

guidata(hObject,data);