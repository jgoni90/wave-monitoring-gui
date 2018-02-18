function UserFcn_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if ~isempty(data.plotid)
    
    %Create User Defined Wave
    ifcn = str2func(hObject.String{hObject.Value});
    UserWave = ifcn(data);
    
    %Update new wave
    set(data.plotid(1),'FaceVertexCData',UserWave)
    data.plotid(1).Vertices(:,3) = UserWave;
    
    %Update maximum, minimum and average values
    UpdateValues(handles);
    
    %Delete pointers, if they exist
    if ~isempty(data.Pointer)
        delete(data.Pointer)
        data.Pointer = [];
    end
    
    %Update visualization options
%     set(handles.View_Plot,'Value',1)
%     set(handles.View_PlotTransparency,'Value',1)
    
end
end