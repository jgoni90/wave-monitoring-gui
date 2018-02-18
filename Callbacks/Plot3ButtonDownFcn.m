function Plot3ButtonDownFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if data.ButtonFunctionsFlag == true;
    data.ButtonDownPoint = eventdata.IntersectionPoint;
    %Delete previous plots, if they exist
    delete(data.Pointer)
    %Draw pointer
    data.Pointer = Draw3DPoint(handles.axesPlot3,data.ButtonDownPoint,data.plotid(3),'IP','{\rmm}^{-1}');
    if eventdata.Button == 3;
        delete(data.Pointer)
        data.Pointer = [];
    end
end

guidata(handles.MainFigure,data);