function View_MaxMin_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = get(hObject,'Value');
if value
    data.Maximum{1} = DrawMax(handles.axesPlot1,data.plotid(1),'U','m');
    data.Minimum{1} = DrawMin(handles.axesPlot1,data.plotid(1),'U','m');
    data.Maximum{2} = DrawMax(handles.axesPlot2,data.plotid(2),'H','m');
    data.Minimum{2} = DrawMin(handles.axesPlot2,data.plotid(2),'H','m');
    data.Maximum{3} = DrawMax(handles.axesPlot3,data.plotid(3),'IP','{\rmm}^{-1}');
    data.Minimum{3} = DrawMin(handles.axesPlot3,data.plotid(3),'IP','{\rmm}^{-1}');
    data.Maximum{4} = DrawMax(handles.axesPlot4,data.plotid(4),'Pot.','W/m');
    data.Minimum{4} = DrawMin(handles.axesPlot4,data.plotid(4),'Pot.','W/m');
else
    delete(data.Maximum)
    delete(data.Minimum)
end

guidata(handles.MainFigure,data);