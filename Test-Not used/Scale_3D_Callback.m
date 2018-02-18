function Scale_3D_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

% New scale value
value = str2double(get(hObject,'String'));
if isnan(value)
    value = 1;
end
set(hObject,'String',num2str(value))

% %Draw solution
% if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
%     if ~isempty(data.plotid(1))
%         cla(handles.axesPlot1)
%         delete(data.plotid(1))
%         delete(data.bottomid(1))
%     end
%     [data.plotid(1),data.bottomid(1)] = UpdatePlot1(handles,value);
% elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
%     handles.axesPlot2.DataAspectRatio(3) = 1/value;
% end
DoPlot(data,handles,false);

guidata(handles.MainFigure,data);