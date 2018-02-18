function View_MeshTransparency_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = hObject.Value;

set(data.meshid,'FaceAlpha',value)

end