function Info_Panel_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if ~isempty(data.rectSelectionid)
    limits = [min(get(data.rectSelectionid,'XData')) max(get(data.rectSelectionid,'XData'))...
        min(get(data.rectSelectionid,'YData')) max(get(data.rectSelectionid,'YData'))];
    pos = createRuledConnectivity(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,limits);
    Trect = data.PGDmeshes(1).T.int(pos,:);
    if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
        axes = handles.axesPlot1;
        plot = data.plotid(1);
    elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
        axes = handles.axesPlot2;
        plot = data.plotid(2);
    elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3')
        axes = handles.axesPlot3;
        plot = data.plotid(3);
    elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4')
        axes = handles.axesPlot4;
        plot = data.plotid(4);
    end
    if get(handles.Average_Text,'Value') == 1;
        candidates = plot.Vertices(Trect,:);
        idx = knnsearch(candidates(:,3),mean(plot.Vertices(Trect,3)));
        xy = candidates(idx,1:2);
        if ~isempty(data.Pointer)
            set(data.Pointer(1),'XData',xy(1),'YData',xy(2))
        end
        [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
    elseif get(handles.Minimum_Text,'Value') == 1;
        candidates = plot.Vertices(Trect,:);
        idx = knnsearch(candidates(:,3),min(plot.Vertices(Trect,3)));
        xy = candidates(idx,1:2);
        if ~isempty(data.Pointer)
            set(data.Pointer(1),'XData',xy(1),'YData',xy(2))
        end
        [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
    elseif get(handles.Maximum_Text,'Value') == 1;
        candidates = plot.Vertices(Trect,:);
        idx = knnsearch(candidates(:,3),max(plot.Vertices(Trect,3)));
        xy = candidates(idx,1:2);
        if ~isempty(data.Pointer)
            set(data.Pointer(1),'XData',xy(1),'YData',xy(2))
        end
        [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
    end
end

guidata(handles.MainFigure,data);