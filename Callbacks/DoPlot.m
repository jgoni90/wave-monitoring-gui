function [plots,bottoms,background,mesh,boundaries,colorbars] = DoPlot(data,handles,update)

% General solution
% Frequencies = data.PGDmeshes(2).X; % 50x1 array
% Directions = data.PGDmeshes(3).X; %50x1 array
if ~update
%     U = interpolatePGD(data.pgd,2,{[data.PGDmeshes(2).X,data.PGDmeshes(3).X]},data.OPTmeshes,data.nOfPGDterms);
    ifcn = str2func(handles.UserFcn.String{handles.UserFcn.Value});
    UserWave = ifcn(data); %Definition of generalized wave
end

%Interpolate PGD
% upgd = interpolatePGD(data.pgd,2:data.parameters.nOfPGDdimensions,...
%     data.snapshot,data.PGDmeshes,data.nOfPGDterms);
upgd = interpolatePGD(data.pgd,2,...
    {[data.snapshot{1},data.snapshot{2}]},data.OPTmeshes,data.nOfPGDterms);

%Protective index
IP = 1./(abs(upgd));

Pot = (9.8^2*9800.*abs(upgd).^2*(2*pi)/data.snapshot{1})/(32*pi); %Definition of potential

%Plot PGD solution
if update
    set(data.plotid(2),'FaceVertexCData',abs(upgd))
    data.plotid(2).Vertices(:,3) = abs(upgd);
    set(data.plotid(3),'FaceVertexCData',IP)
    data.plotid(3).Vertices(:,3) = IP;
    set(data.plotid(4),'FaceVertexCData',Pot)
    data.plotid(4).Vertices(:,3) = Pot;
else
    
    % Bottom values
    if ~isempty(data.bottom)
        bottom = data.bottom;
    else
        bottom = -15.0119*ones(size(data.PGDmeshes(1).X,1),1);
    end
    
    axes(handles.axesPlot1) %Tab 1
    axis tight
    background(1) = imagesc(data.background);
    set(handles.axesPlot1,'YDir','normal','Visible','off')
	plots(1) = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        UserWave,data.PGDmeshes(1).referenceElement);
    caxis auto
    [bottoms(1),colorbars(1)] = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        bottom,data.PGDmeshes(1).referenceElement);
    bottoms(1).Visible = 'off';
    mesh(1) = plotMesh(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,data.PGDmeshes(1).referenceElement.faceNodes);
    mesh(1).Visible = 'off';
    boundaries{1} = plotBoundaries(data);
%     XLim = [min(min(bottoms(1).XData)) max(max(bottoms(1).XData))];
%     YLim = [min(min(bottoms(1).YData)) max(max(bottoms(1).YData))];
    XLim = [min(data.PGDmeshes(1).X(:,1)) max(data.PGDmeshes(1).X(:,1))];
    YLim = [min(data.PGDmeshes(1).X(:,2)) max(data.PGDmeshes(1).X(:,2))];
    set(background,'XData',XLim)
    set(background,'YData',YLim)
    
    axes(handles.axesPlot2) %Tab 2
    background(2) = imagesc(data.background);
    set(handles.axesPlot2,'YDir','normal','Visible','off')
	plots(2) = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        abs(upgd),data.PGDmeshes(1).referenceElement);
    caxis auto
    [bottoms(2),colorbars(2)] = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        bottom,data.PGDmeshes(1).referenceElement);
    bottoms(2).Visible = 'off';
    mesh(2) = plotMesh(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,data.PGDmeshes(1).referenceElement.faceNodes);
    mesh(2).Visible = 'off';
    boundaries{2} = plotBoundaries(data);
    set(background,'XData',XLim)
    set(background,'YData',YLim)

    axes(handles.axesPlot3) %Tab 3
    background(3) = imagesc(data.background);
    set(handles.axesPlot3,'YDir','normal','Visible','off')
	plots(3) = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        IP,data.PGDmeshes(1).referenceElement);
    caxis auto
    [bottoms(3),colorbars(3)] = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        bottom,data.PGDmeshes(1).referenceElement);
    bottoms(3).Visible = 'off';
    mesh(3) = plotMesh(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,data.PGDmeshes(1).referenceElement.faceNodes);
    mesh(3).Visible = 'off';
    boundaries{3} = plotBoundaries(data);
    set(background,'XData',XLim)
    set(background,'YData',YLim)
    
    axes(handles.axesPlot4) % Tab 4
    background(4) = imagesc(data.background);
    set(handles.axesPlot4,'YDir','normal','Visible','off')
	plots(4) = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        Pot,data.PGDmeshes(1).referenceElement);
    caxis auto
    [bottoms(4),colorbars(4)] = plot3DSolution(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,...
        bottom,data.PGDmeshes(1).referenceElement);
    bottoms(4).Visible = 'off';
    mesh(4) = plotMesh(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,data.PGDmeshes(1).referenceElement.faceNodes);
    mesh(4).Visible = 'off';
    boundaries{4} = plotBoundaries(data);
    set(background,'XData',XLim)
    set(background,'YData',YLim)
        
    axes(handles.axesArrow) %Change current axes
end