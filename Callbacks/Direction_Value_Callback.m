function Direction_Value_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

%New frequency value
value = str2double(get(hObject,'String'));
maxdir = 180*data.PGDmeshes(3).X(end)/pi;
mindir = 180*data.PGDmeshes(3).X(1)/pi;
if value > maxdir
    value = maxdir;
elseif value < mindir
    value = mindir;
end
set(hObject,'String',num2str(value))
data.snapshot{2} = pi*value/180;
set(data.lineid{2},'XData',[180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2}])
set(data.rectid{2},'YData',[180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2}])
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
data.arrowdata.p0p(1) = data.snapshot{2} - pi;
[data.arrowdata.p0(1),data.arrowdata.p0(2)] = ...
    pol2cart(data.arrowdata.p0p(1),data.arrowdata.p0p(2));

%Draw arrow and solution
delete(data.arrowid)
data.arrowid = DoArrow(data,handles);
DoPlot(data,handles,true);

guidata(handles.MainFigure,data);