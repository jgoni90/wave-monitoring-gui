function Plot2ButtonDownFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if data.ButtonFunctionsFlag == true;
    data.ButtonDownPoint = eventdata.IntersectionPoint;
    %Delete previous plots, if they exist
    delete(data.Pointer)
    %Draw pointer
    data.Pointer = Draw3DPoint(handles.axesPlot2,data.ButtonDownPoint,data.plotid(2),'H','m');    
    if eventdata.Button == 3;
        delete(data.Pointer)
        data.Pointer = [];
    end
end

guidata(handles.MainFigure,data);