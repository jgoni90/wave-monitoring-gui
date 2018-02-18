function View_PlotTransparency_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = hObject.Value;

set(data.plotid,'FaceAlpha',value)

end