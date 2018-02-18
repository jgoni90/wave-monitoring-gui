function PanXY_Callback(hObject,eventdata,handles)

switch get(hObject,'State')
    case 'on'
        set([handles.ZoomXY handles.Rotate],'State','off')
        pan(handles.axesPlot1,'on')
    case 'off'
        pan(handles.axesPlot1,'off')
end