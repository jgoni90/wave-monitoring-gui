function Colorbar2ButtonDownFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if strcmp(get(handles.MainFigure,'SelectionType'),'normal')
    limit = eventdata.IntersectionPoint(2);

    if ~isempty(data.limitid)
        delete(data.limitid)
        data.limitid = [];
    end
    if ~isempty(data.limitlineid)
        delete(data.limitlineid{1})
        delete(data.limitlineid{2})
        data.limitlineid = [];
    end
    if ~isempty(data.limitrectid)
        delete(data.limitrectid)
        data.limitrectid = [];
    end

    axes(handles.axesPlot2)

    % contourHandle = contour(data.plotid(1).Vertices,[limit limit]);

    data.limitid = patch('XData',[handles.axesPlot2.XLim(1) handles.axesPlot2.XLim(2) handles.axesPlot2.XLim(2) handles.axesPlot2.XLim(1)],...
        'YData',[handles.axesPlot2.YLim(1) handles.axesPlot2.YLim(1) handles.axesPlot2.YLim(2) handles.axesPlot2.YLim(2)],...
        'ZData',[limit limit limit limit],...
        'FaceColor','k','FaceAlpha',0.25,'EdgeAlpha',0,'Parent',handles.axesPlot2);

    if ~isempty(data.graphid)
        hold on
        data.limitlineid{1} = line(handles.axesTheta2D.XLim,[limit limit],...
            'Color','k','Parent',handles.axesTheta2D);
        data.limitlineid{2} = line(handles.axesOmega2D.XLim,[limit limit],...
            'Color','k','Parent',handles.axesOmega2D);
        hold off
    end

    if ~isempty(data.surfid)
        hold on
        data.limitrectid = patch('XData',[handles.axes3D.XLim(1) handles.axes3D.XLim(2) handles.axes3D.XLim(2) handles.axes3D.XLim(1)],...
            'YData',[handles.axes3D.YLim(1) handles.axes3D.YLim(1) handles.axes3D.YLim(2) handles.axes3D.YLim(2)],...
            'ZData',[limit limit limit limit],'FaceColor','k','FaceAlpha',0.25,'EdgeAlpha',0,'Parent',handles.axes3D);
        hold off
    end
elseif strcmp(get(handles.MainFigure,'SelectionType'),'open')
    if ~isempty(data.limitid)
        delete(data.limitid)
        data.limitid = [];
    end
    if ~isempty(data.limitlineid)
        delete(data.limitlineid{1})
        delete(data.limitlineid{2})
        data.limitlineid = [];
    end
    if ~isempty(data.limitrectid)
        delete(data.limitrectid)
        data.limitrectid = [];
    end
end
    
guidata(handles.MainFigure,data)