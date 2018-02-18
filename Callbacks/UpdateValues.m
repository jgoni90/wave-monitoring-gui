function UpdateValues(handles)

data = guidata(handles.MainFigure);

precision = 3;

if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
    set(handles.Average_Text,'String','Average U')
    set(handles.Average_Value,'String',[num2str(mean(data.plotid(1).Vertices(:,3)),precision) ' m'])
    set(handles.Minimum_Text,'String','Minimum U')
    set(handles.Minimum_Value,'String',[num2str(min(data.plotid(1).Vertices(:,3)),'%.3e') ' m'])
    set(handles.Maximum_Text,'String','Maximum U')
    set(handles.Maximum_Value,'String',[num2str(max(data.plotid(1).Vertices(:,3)),precision) ' m'])
elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
    set(handles.Average_Text,'String','Average H')
    set(handles.Average_Value,'String',[num2str(mean(data.plotid(2).Vertices(:,3)),precision) ' m'])
    set(handles.Minimum_Text,'String','Minimum H')
    set(handles.Minimum_Value,'String',[num2str(min(data.plotid(2).Vertices(:,3)),'%.3e') ' m'])
    set(handles.Maximum_Text,'String','Maximum H')
    set(handles.Maximum_Value,'String',[num2str(max(data.plotid(2).Vertices(:,3)),precision) ' m'])
elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3')
    set(handles.Average_Text,'String','Average IP')
    set(handles.Average_Value,'String',[num2str(mean(data.plotid(3).Vertices(:,3)),precision) ' [1/m]'])
    set(handles.Minimum_Text,'String','Minimum IP')
    set(handles.Minimum_Value,'String',[num2str(min(data.plotid(3).Vertices(:,3)),precision) ' [1/m]'])
    set(handles.Maximum_Text,'String','Maximum IP')
    set(handles.Maximum_Value,'String',[num2str(max(data.plotid(3).Vertices(:,3)),precision) ' [1/m]'])
elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4')
    set(handles.Average_Text,'String','Average P')
    set(handles.Average_Value,'String',[num2str(mean(data.plotid(4).Vertices(:,3))/1000,precision) ' kW/m'])
    set(handles.Minimum_Text,'String','Minimum P')
    set(handles.Minimum_Value,'String',[num2str(min(data.plotid(4).Vertices(:,3))/1000,'%.3e') ' kW/m'])
    set(handles.Maximum_Text,'String','Maximum P')
    set(handles.Maximum_Value,'String',[num2str(max(data.plotid(4).Vertices(:,3))/1000,precision) ' kW/m'])
end

guidata(handles.MainFigure,data);