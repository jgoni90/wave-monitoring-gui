function MainFigureWindowButtonDownFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if data.ButtonFunctionsFlag == true;

%% Determine current axes
PlotAxesPos = handles.MainTabs.Position;
ArrowAxesPos = [handles.Input_Panel.Position(1)+handles.Input_Panel.Position(3)*handles.axesArrow.Position(1)...
    handles.Input_Panel.Position(2)+handles.Input_Panel.Position(4)*handles.axesArrow.Position(2)...
    handles.Input_Panel.Position(3)*handles.axesArrow.Position(3)...
    handles.Input_Panel.Position(4)*handles.axesArrow.Position(4)];
CurrentPoint = get(handles.MainFigure,'CurrentPoint');

ArrowX = ArrowAxesPos(1) <= CurrentPoint(1) && CurrentPoint(1) <= ArrowAxesPos(1) + ArrowAxesPos(3);
ArrowY = ArrowAxesPos(2) <= CurrentPoint(2) && CurrentPoint(2) <= ArrowAxesPos(2) + ArrowAxesPos(4);

PlotX = PlotAxesPos(1) <= CurrentPoint(1) && CurrentPoint(1) <= PlotAxesPos(1) + PlotAxesPos(3);
PlotY = PlotAxesPos(2) <= CurrentPoint(2) && CurrentPoint(2) <= PlotAxesPos(2) + PlotAxesPos(4);

if ArrowX && ArrowY
    %Current point in the figure
    cp = get(handles.axesArrow,'CurrentPoint');
    xy = [cp(1,1) cp(1,2)];

    %Selection point
    p0 = data.arrowdata.p0;
    circle = (xy(1)-p0(1))^2 + (xy(2)-p0(2))^2;
    circlerad = pi/15;

    %Set properties
    if circle <= circlerad^2
        data.userSelection = true;
        set(data.arrowid(3),'markerfacecolor','r','markeredgecolor','r')
    end
elseif PlotX && PlotY
    if strcmp(get(handles.MainFigure,'SelectionType'),'open')
        data.rectSelection = false;
        %Delete previous plots, if they exist
        if ~isempty(data.graphid)
            delete(data.graphid{1});
            delete(data.graphid{2});
            data.graphid = [];
        end
        xy = data.ButtonDownPoint(1:2);
        [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
    elseif strcmp(get(handles.MainFigure,'SelectionType'),'normal')
        if ~isempty(data.rectSelectionid)
            delete(data.rectSelectionid);
            data.rectSelectionid = [];
        end
        data.rectSelection = true;
        data.rectIniPoint = [handles.axesPlot1.CurrentPoint(1,1) handles.axesPlot1.CurrentPoint(1,2)];
        UpdateValues(handles);
    end
end
end

guidata(handles.MainFigure,data);