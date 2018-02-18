function View_BGTransparency_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = hObject.Value;

set(data.backgroundid,'AlphaData',value)

end