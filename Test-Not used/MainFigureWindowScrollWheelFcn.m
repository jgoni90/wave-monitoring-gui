function MainFigureWindowScrollWheelFcn(hObject,eventdata,handles)

%Currently not fully implemented (needs tweaking,use ZoomXY_Callback
%instead)

% data = guidata(handles.MainFigure);
% 
% if data.ButtonFunctionsFlag == true;
% 
% axes(handles.axesPlot1)
% 
% %% Determine current axes
% PlotAxesPos = get(handles.MainTabs,'Position');
% CurrentPoint = get(handles.MainFigure,'CurrentPoint');
% 
% PlotX = PlotAxesPos(1) <= CurrentPoint(1) && CurrentPoint(1) <= PlotAxesPos(1) + PlotAxesPos(3);
% PlotY = PlotAxesPos(2) <= CurrentPoint(2) && CurrentPoint(2) <= PlotAxesPos(2) + PlotAxesPos(4);
% 
% XLim = get(handles.axesPlot1,'XLim');
% YLim = get(handles.axesPlot1,'YLim');
% PlotDX = diff(get(handles.axesPlot1,'XLim'));
% PlotDY = diff(get(handles.axesPlot1,'YLim'));
% Cx =  XLim(1) + (CurrentPoint(1)-PlotAxesPos(1)).*(PlotDX/PlotAxesPos(3));
% Cy =  YLim(1) + (CurrentPoint(2)-PlotAxesPos(2)).*(PlotDY/PlotAxesPos(4));
% 
% axlims = axis;
% zf = 1/8;     %increment of zoom per mouse scroll-wheel click
% axsize = [abs(axlims(2)-axlims(1)),abs(axlims(4)-axlims(3))];
% 
% if PlotX && PlotY
%     if eventdata.VerticalScrollCount > 0           %Zoom out
%       if ~isempty(CurrentPoint)
%           axlims = [Cx-(axsize(1)*(1+zf)/2),...    %xmin
%               Cx+(axsize(1)*(1+zf)/2),...
%               Cy-(axsize(2)*(1+zf)/2),...
%               Cy+(axsize(2)*(1+zf)/2)];
%       else axlims = [axlims(1)-axsize(1)*zf,axlims(2)+axsize(1)*zf,...
%               axlims(3)-axsize(2)*zf,axlims(4)+axsize(2)*zf];
%       end
% 
%     elseif eventdata.VerticalScrollCount < 0       %Zoom in
%       if ~isempty(CurrentPoint)
%           axlims = [Cx-(axsize(1)*(1-zf)/2),...    %xmin
%               Cx+(axsize(1)*(1-zf)/2),...
%               Cy-(axsize(2)*(1-zf)/2),...
%               Cy+(axsize(2)*(1-zf)/2)];
%       else axlims = [axlims(1)+axsize(1)*zf,axlims(2)-axsize(1)*zf,...
%               axlims(3)+axsize(2)*zf,axlims(4)-axsize(2)*zf];
%       end
%     end
%     
% if abs(data.PlotAxesCoordinates) < abs(axlims);
%     axis(data.PlotAxesCoordinates);
% else
%     axis(axlims);
% end
%     
% end
% end

% guidata(handles.MainFigure,data);