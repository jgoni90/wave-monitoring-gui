function View_BotTransparency_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = hObject.Value;

set(data.bottomid,'FaceAlpha',value)

end