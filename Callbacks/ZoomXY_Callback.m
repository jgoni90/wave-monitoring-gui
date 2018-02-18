function ZoomXY_Callback(hObject,eventdata,handles)

switch get(hObject,'State')
    case 'on'
        set([handles.PanXY handles.Rotate],'State','off')
        zoom(handles.axesPlot1,'on')
        zoomHandle = zoom(handles.axesPlot1);
        set(zoomHandle,'RightClickAction','InverseZoom')
    case 'off'
        zoom(handles.axesPlot1,'off')
end
        