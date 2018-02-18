function MainFigureWindowButtonMotionFcn(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

if data.ButtonFunctionsFlag == true;
    
    %% Determine current axes
    PlotAxesPos = get(handles.MainTabs,'Position');
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
        if data.userSelection
            if ~isempty(data.rectSelectionid)
                delete(data.rectSelectionid)
                data.rectSelectionid = [];
            end
            
            %Current point in the figure
            cp = get(handles.axesArrow,'currentPoint');
            data.arrowdata.p0 = [cp(1,1) cp(1,2)];
            
            %Take snapshot and corrected point
            [p0,p0p,freq,theta] = getParamValuesFromCP(data);
            data.arrowdata.prev_p0 = data.arrowdata.p0;
            data.arrowdata.prev_p0p = data.arrowdata.p0p;
            data.arrowdata.p0 = p0;
            data.arrowdata.p0p = p0p;
            data.snapshot{1} = freq;
            data.snapshot{2} = theta + pi;
            
            %Draw
            try delete(data.arrowid), catch return, end
            data.arrowid = DoArrow(data,handles);
            %     DoPlot(data,handles,true);
            
            %Update selection panel
            set(handles.Freq_Value,'String',num2str(data.snapshot{1}))
            set(handles.Period_Value,'String',num2str(2*pi/data.snapshot{1}))
            set(handles.Direction_Value,'String',num2str(180*data.snapshot{2}/pi))
            
            %Update 2D plot lines
            if ~isempty(data.lineid)
                set(data.lineid{1},'XData',[2*pi/data.snapshot{1} 2*pi/data.snapshot{1}])
                set(data.rectid{1},'XData',[2*pi/data.snapshot{1} 2*pi/data.snapshot{1} 2*pi/data.snapshot{1} 2*pi/data.snapshot{1}])
                set(data.lineid{2},'XData',[180*data.snapshot{2}/pi 180*data.snapshot{2}/pi])
                set(data.rectid{2},'YData',[180*data.snapshot{2}/pi 180*data.snapshot{2}/pi 180*data.snapshot{2}/pi 180*data.snapshot{2}/pi])
            end
            
            %Update Pointer
            if size(data.Pointer,2) ~= 1
                delete(data.Pointer(2))
                data.Pointer(2) = [];
            end
        end
    elseif PlotX && PlotY
        if strcmp(get(handles.XCoord_Value,'Enable'),'on') && strcmp(get(handles.YCoord_Value,'Enable'),'on') 
            handles.XCoord_Value.String = handles.axesPlot1.CurrentPoint(1,1);
            handles.YCoord_Value.String = handles.axesPlot1.CurrentPoint(1,2);
            if ~isempty(data.rectSelection) && data.rectSelection
                if ~isempty(data.Pointer)
                    delete(data.Pointer)
                    data.Pointer = [];
                end
                % get the top left corner and width and height of 
                % rectangle (note the absolute value forces it to "open"
                % from left to right - need smarter logic for other direction)
                xi = data.rectIniPoint(1);
                yi = data.rectIniPoint(2);
                xf = handles.axesPlot1.CurrentPoint(1,1);
                yf = handles.axesPlot1.CurrentPoint(1,2);
                if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1')
                    z = max(data.plotid(1).Vertices(:,3));
                    parent = handles.axesPlot1;
                elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2')
                    z = max(data.plotid(2).Vertices(:,3));
                    parent = handles.axesPlot2;
                elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3')
                    z = max(data.plotid(3).Vertices(:,3));
                    parent = handles.axesPlot3;
                elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4')
                    z = max(data.plotid(4).Vertices(:,3));
                    parent = handles.axesPlot4;
                end
                % rectangle drawn in white (better colour needed for different
                % images?)
                if isempty(data.rectSelectionid)
                    % empty so rectangle not yet drawn
                    data.rectSelectionid = patch('XData',[xi xf xf xi],'YData',[yi yi yf yf],'ZData',[z z z z],...
                        'EdgeColor','r','LineStyle','--','FaceColor','r','FaceAlpha',0.5,'Parent',parent);
                else
                    % need to redraw
                    set(data.rectSelectionid,'XData',[xi xf xf xi],'YData',[yi yi yf yf]);
                end
            end
        end
    end
end

guidata(handles.MainFigure,data);