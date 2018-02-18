function load_pgd(hObject,eventdata,handles)

data = guidata(hObject);
directory = [pwd '/PGD_Data'];

set(data.PortMarkers,'Color','k','MarkerFaceColor','k')
set(hObject,'Color','r','MarkerFaceColor','r')
delete(data.PortInfo)

if strcmp(get(handles.MainFigure,'SelectionType'),'open')
    
    %Select firs tab as current
    set(handles.MainTabs,'SelectedTab',handles.Tab1)
    
    %Reset axes
    cla(handles.axesArrow)
    cla(handles.axesPlot1)
    cla(handles.axesPlot2)
    cla(handles.axesPlot3)
    cla(handles.axesPlot4)
    cla(handles.axes3D)
    cla(handles.axesOmega2D)
    cla(handles.axesTheta2D)
    
    %Delete previous plots, if they exist
%     if ~isempty(data.plotid)
%         delete(data.plotid)
%         data.plotid = [];
%     end
%     if ~isempty(data.arrowid)
%         delete(data.arrowid)
%         data.arrowid = [];
%     end
%     if ~isempty(data.backgroundid)
%         delete(data.backgroundid)
%         data.backgroundid = [];
%     end
%     if ~isempty(data.bottomid)
%         delete(data.bottomid)
%         data.bottomid = [];
%     end
%     if ~isempty(data.meshid)
%         delete(data.meshid)
%         data.meshid = [];
%     end
%     if ~isempty(data.boundaryid)
%         delete(data.boundaryid)
%         data.boundaryid = [];
%     end
%     if ~isempty(data.Pointer)
%         delete(data.Pointer)
%         data.Pointer = [];
%     end
%     if ~isempty(data.graphid)
%         delete(data.graphid)
%         data.graphid = [];
%     end
%     if ~isempty(data.colorbarid)
%         delete(data.colorbarid)
%         data.colorbarid = [];
%     end
%     if ~isempty(data.lineid)
%         delete(data.lineid)
%         data.lineid = [];
%     end
%     if ~isempty(data.surfid)
%         delete(data.surfid)
%         data.surfid = [];
%     end
%     if ~isempty(data.rectid)
%         delete(data.rectid)
%         data.rectid = [];
%     end
%     if ~isempty(data.limitid)
%         delete(data.limitid)
%         data.limitid = [];
%     end
    
    %Locate files to load
    dir_path = [directory '/' get(hObject,'Tag')];
    bottom_dir_path = [dir_path '/Bottom'];
    image = dir([dir_path '/*.png']);
    files = dir([dir_path '/*.mat']);
    bottom = dir([bottom_dir_path '/*.mat']);
    filelist = {files.name}';
    
    % Load PGD Data
    setOutput({'Loading PGD data...'},handles.Text_Output,1)
    for i = 1:size(filelist,1)
        load(filelist{i})
    end
    data.background = imread(image.name);
    bottomfile = bottom.name;
    load(bottomfile);
    data.bottom = bottom;
    
        data.pgd           = projectedpgd;
        data.PGDalgorithm  = PGDalgorithm;
        data.parameters    = parameters;
        data.PGDmeshes     = PGDmeshes;
        data.OPTmeshes(2)  = PROJmeshes;
        data.parameters.nOfPGDdimensions = data.parameters.nOfPGDdimensions - 1;
        
        data.solution = data.pgd.RB{1}*data.pgd.RB{2}';
        
        %Initial parameters configuration (snapshot & #PGD terms)
        
        %Arrow fixed data
        fac = PGDmeshes(2).X(end) - PGDmeshes(2).X(1);
        data.arrowdata.a0 = fac / pi;
        data.arrowdata.a1 = 13 / fac;
        data.arrowdata.b0 = PGDmeshes(2).X(1);
        data.arrowdata.b1 = 2 - 13*PGDmeshes(2).X(1) / fac;
        data.arrowdata.thetalim0 = data.PGDmeshes(3).X(1) - pi;
        data.arrowdata.thetalim1 = data.PGDmeshes(3).X(end) - pi;
        
        %Arrow variable data
        data.arrowdata.p0 = [0,0];
        data.arrowdata.prev_p0 = [0,0];
        data.arrowdata.p0p = [0,0];
        data.arrowdata.prev_p0p = [0,0];
        
        %Initial snapshot
        data.snapshot = cell(data.parameters.nOfPGDdimensions,1);
        
            %Freq and period
            data.snapshot{1} = data.PGDmeshes(2).X(1);
            set(handles.Freq_Value,'String',num2str(data.snapshot{1}))
            set(handles.Period_Value,'String',num2str(2*pi/data.snapshot{1}))

            %Incoming wave Direction
            data.snapshot{2} = data.PGDmeshes(3).X(1);
            set(handles.Direction_Value,'String',num2str(180*data.snapshot{2}/pi))

            %Rest of dimensions (ALPHA dims)
            for idim = 4:data.parameters.nOfPGDdimensions
                data.snapshot{idim} = data.PGDmeshes(idim+1).X(1);
            end

            %Number of PGD terms
            data.nOfPGDterms = size(data.pgd.RB{1},2);
            set(handles.PGDTerms_Value,'String',num2str(data.nOfPGDterms))
            
            %Element degree
            set(handles.ElemType_Value,'String',sprintf('Order %d Tri',...
                (data.PGDmeshes(1).referenceElement.degree)))
            
            %Periods
            set(handles.PeriodRange_Value,'String',sprintf('%d-%d s',...
                (data.parameters.meshes(2).Tini),(data.parameters.meshes(2).Tend)))
            
            %Directions
            set(handles.DirectionRange_Value,'String',sprintf('%d-%d%c',...
                (data.parameters.meshes(3).THETAini*180/pi()),...
                (data.parameters.meshes(3).THETAend*180/pi()),char(176))) %Minimum direction
        
        %Do initial plots
        setOutput({'Plotting figure and arrow data...'},handles.Text_Output)
        [data.plotid,data.bottomid,data.backgroundid,data.meshid,...
            data.boundaryid,data.colorbarid] = DoPlot(data,handles,false);
        data.arrowid = DoArrow(data,handles);
        
        %Functions
        data.plotid(1).ButtonDownFcn = {@Plot1ButtonDownFcn,handles};
        data.plotid(2).ButtonDownFcn = {@Plot2ButtonDownFcn,handles};
        data.plotid(3).ButtonDownFcn = {@Plot3ButtonDownFcn,handles};
        data.plotid(4).ButtonDownFcn = {@Plot4ButtonDownFcn,handles};
        
        %Functions for colorbarhandles
        set(data.colorbarid(1),'ButtonDownFcn',{@Colorbar1ButtonDownFcn,handles})
        data.colorbarid(1).UIContextMenu.Visible = 'off';
        data.colorbarid(2).ButtonDownFcn = {@Colorbar2ButtonDownFcn,handles};
        data.colorbarid(3).ButtonDownFcn = {@Colorbar3ButtonDownFcn,handles};
        data.colorbarid(4).ButtonDownFcn = {@Colorbar4ButtonDownFcn,handles};
        
        %Load boundary list
        BCList = fieldnames(data.PGDmeshes(1).elementFaceInfo);
        set(handles.View_BC,'String',{'None',BCList{:},'All'})
        
        %Maximum, minimum and average
        set(handles.Average_Text,'String','Average U')
        set(handles.Average_Value,'String',[num2str(mean(data.plotid(1).Vertices(:,3)),3) ' m'])
        set(handles.Minimum_Text,'String','Minimum U')
        set(handles.Minimum_Value,'String',[num2str(min(data.plotid(1).Vertices(:,3)),'%.3e') ' m'])
        set(handles.Maximum_Text,'String','Maximum U')
        set(handles.Maximum_Value,'String',[num2str(max(data.plotid(1).Vertices(:,3)),3) ' m'])
         
        set([handles.ZoomXY handles.PanXY handles.Rotate handles.View_Mesh...
            handles.View_Bottom handles.View_Plot handles.View_Background...
            handles.View_BC handles.XCoord_Value handles.YCoord_Value handles.View_BC...
            handles.View_PlotTransparency handles.View_BGTransparency...
            handles.Average_Value handles.Minimum_Value handles.Maximum_Value...
            handles.Average_Text handles.Minimum_Text handles.Maximum_Text],'Enable','on')
        set([handles.View_Background handles.View_Plot handles.View_BC handles.View_BotTransparency...
            handles.View_BGTransparency handles.View_MeshTransparency handles.View_PlotTransparency...
            handles.Average_Text],'value',1)
        set([handles.View_Axes handles.View_Bottom handles.View_Mesh],'value',0)
        % Synchronize axes' limits
        uistack(handles.Input_Panel,'bottom')
        data.ButtonFunctionsFlag = true;
        data.PlotAxesCoordinates = axis(handles.axesPlot1);
        guidata(handles.MainFigure,data);
        setOutput({'Done!'},handles.Text_Output)
elseif strcmp(get(handles.MainFigure,'SelectionType'),'normal')
    index_selected = strmatch(get(hObject,'Tag'),data.PortNames);
    offset = 25;
    data.PortInfo = text(handles.axesLocation.XLim(1)+offset,handles.axesLocation.YLim(2)-70,...
        {[get(hObject,'Tag')],['Lat. ',num2str(data.PortCoordinates(index_selected,1))],...
        ['Long. ',num2str(data.PortCoordinates(index_selected,2))]},'BackgroundColor','w','Parent',handles.axesLocation);
    guidata(handles.MainFigure,data);
end

end