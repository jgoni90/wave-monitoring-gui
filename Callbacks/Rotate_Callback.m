function Rotate_Callback(hObject,eventdata,handles)

% linkprop([handles.axesPlot1 handles.axesPlot2 handles.axesPlot3 handles.axesPlot4],'CameraPosition');

switch get(hObject,'State')
    case 'on'
        set([handles.PanXY handles.ZoomXY],'State','off')
        if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
            rotate3d(handles.axesPlot1,'on')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
            rotate3d(handles.axesPlot2,'on')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3')
            rotate3d(handles.axesPlot3,'on')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4')
            rotate3d(handles.axesPlot4,'on')
        end
    case 'off'
        if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
            rotate3d(handles.axesPlot1,'off')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
            rotate3d(handles.axesPlot2,'off')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3')
            rotate3d(handles.axesPlot3,'off')
        elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4')
            rotate3d(handles.axesPlot4,'off')
        end
end