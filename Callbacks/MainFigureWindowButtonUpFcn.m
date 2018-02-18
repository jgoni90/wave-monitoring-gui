function MainFigureWindowButtonUpFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if data.userSelection
    data.userSelection = false;
    set(data.arrowid(3),'markerfacecolor','k','markeredgecolor','k')

    %Draw solution
    DoPlot(data,handles,true);
    UpdateValues(handles);
elseif data.rectSelection
    % reset the button down flag
    data.rectSelection = false;
    if ~isempty(data.rectSelectionid)
        set(data.rectSelectionid,'LineStyle','-')
        limits = [min(get(data.rectSelectionid,'XData')) max(get(data.rectSelectionid,'XData'))...
            min(get(data.rectSelectionid,'YData')) max(get(data.rectSelectionid,'YData'))];
        pos = createRuledConnectivity(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,limits);
        Trect = data.PGDmeshes(1).T.int(pos,:);
        if size(Trect,1) ~= 0
            UpdateAreaValues(handles,Trect);
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
                hold on
                data.Pointer(1) = plot3(axes,xy(1),xy(2),max(plot.Vertices(:,3)),'w+',...
                    'Markersize',20,'MarkerFaceColor','w','LineWidth',3,'Visible','on');
                hold off
                [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
            elseif get(handles.Minimum_Text,'Value') == 1;
                candidates = plot.Vertices(Trect,:);
                idx = knnsearch(candidates(:,3),min(plot.Vertices(Trect,3)));
                xy = candidates(idx,1:2);
                hold on
                data.Pointer(1) = plot3(axes,xy(1),xy(2),max(plot.Vertices(:,3)),'w+',...
                    'Markersize',20,'MarkerFaceColor','w','LineWidth',3,'Visible','on');
                hold off
                [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
            elseif get(handles.Maximum_Text,'Value') == 1;
                candidates = plot.Vertices(Trect,:);
                idx = knnsearch(candidates(:,3),max(plot.Vertices(Trect,3)));
                xy = candidates(idx,1:2);
                hold on
                data.Pointer(1) = plot3(axes,xy(1),xy(2),max(plot.Vertices(:,3)),'w+',...
                    'Markersize',20,'MarkerFaceColor','w','LineWidth',3,'Visible','on');
                hold off
                [data.graphid,data.lineid,data.surfid,data.rectid] = DrawPlots(data,handles,xy);
            end
        else
            delete(data.rectSelectionid);
            data.rectSelectionid = [];
        end
    end
end

guidata(handles.MainFigure,data);