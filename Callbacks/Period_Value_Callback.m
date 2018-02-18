function Period_Value_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

%New frequency value
value = str2double(get(hObject,'String'));
maxper = 2*pi/data.PGDmeshes(2).X(1);
minper = 2*pi/data.PGDmeshes(2).X(end);
if value > maxper
    value = maxper;
elseif value < minper
    value = minper;
end
set(hObject,'String',num2str(value))
set(handles.Freq_Value,'String',num2str(2*pi/value))
data.snapshot{1} = 2*pi/value;
set(data.lineid{1},'XData',[2*pi./data.snapshot{1} 2*pi./data.snapshot{1}])
set(data.rectid{1},'XData',[2*pi./data.snapshot{1} 2*pi./data.snapshot{1} 2*pi./data.snapshot{1} 2*pi./data.snapshot{1}])
if ~isempty(data.rectSelectionid)
    delete(data.rectSelectionid)
    data.rectSelectionid = [];
elseif size(data.Pointer,2) ~= 1
    delete(data.Pointer(2))
    data.Pointer(2) = [];
end
UpdateValues(handles);

%Update selection point
data.arrowdata.prev_p0 = data.arrowdata.p0;
data.arrowdata.prev_p0p = data.arrowdata.p0p;
data.arrowdata.p0p(2) = (data.snapshot{1}-data.arrowdata.b0)/data.arrowdata.a0;
[data.arrowdata.p0(1),data.arrowdata.p0(2)] = ...
    pol2cart(data.arrowdata.p0p(1),data.arrowdata.p0p(2));

%Draw arrow and solution
delete(data.arrowid)
data.arrowid = DoArrow(data,handles);
DoPlot(data,handles,true);

guidata(handles.MainFigure,data);