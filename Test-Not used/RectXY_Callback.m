function RectXY_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

switch get(hObject,'State')
    case 'on'
        set(handles.PanXY,'State','off')
        set(handles.ZoomXY,'State','off')
        data.selectionRectangle = getrect(handles.axesPlot1);
        set(handles.RectXY,'State','off')
end

guidata(handles.MainFigure,data);        