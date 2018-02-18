function varargout = DrawPlots(data,handles,xy)

%         e = 50;
%         pos = createRuledConnectivity(data.PGDmeshes(1).X,data.PGDmeshes(1).T.int,[xy(1)-e xy(1)+e xy(2)-e xy(2)+e]);
%         Tpoint = data.PGDmeshes(1).T.int(pos,:);
%         Xpoint = data.PGDmeshes(1).X(Tpoint);
%Presearch not used: Index needed, not position.

%Create auxiliary PGD structure
index = knnsearch(data.PGDmeshes(1).X,xy);
pgdaux = data.pgd;
pgdaux.RB{1} = pgdaux.RB{1}(index,:);
thetaindex = knnsearch(data.OPTmeshes(2).X(:,2),data.snapshot{2},'K',size(data.PGDmeshes(3).X,1));
omegaindex = knnsearch(data.OPTmeshes(2).X(:,1),data.snapshot{1},'K',size(data.PGDmeshes(2).X,1));
pgdauxomega.RB{2} = pgdaux.RB{2}(thetaindex,:);
pgdauxtheta.RB{2} = pgdaux.RB{2}(omegaindex,:);

%Create solution vectors
u_omega = abs(pgdaux.RB{1}*pgdauxomega.RB{2}');
u_theta = abs(pgdaux.RB{1}*pgdauxtheta.RB{2}');
u_omega_theta = abs(pgdaux.RB{1}*data.pgd.RB{2}');
if strcmp(handles.MainTabs.SelectedTab.Tag,'Tab1') || strcmp(handles.MainTabs.SelectedTab.Tag,'Tab2');
    solution_omega = u_omega;
    solution_theta = u_theta;
    solution_omega_theta = u_omega_theta;
elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab3');
    solution_omega = 1./u_omega;
    solution_theta = 1./u_theta;
    solution_omega_theta = 1./u_omega_theta;
elseif strcmp(handles.MainTabs.SelectedTab.Tag,'Tab4');
    solution_omega = (9.8^2*9800.*u_omega.^2*(2*pi)/data.snapshot{1})/(32*pi);
    solution_theta = (9.8^2*9800.*u_theta.^2*(2*pi)/data.snapshot{1})/(32*pi);
    solution_omega_theta = (9.8^2*9800.*u_omega_theta.^2*(2*pi)/data.snapshot{1})/(32*pi);
end
[X,Y] = meshgrid(2*pi./data.PGDmeshes(2).X,180/pi.*data.PGDmeshes(3).X);
x_axis = reshape(X,[size(data.OPTmeshes(2).X,1),1]); %Omega
y_axis = reshape(Y,[size(data.OPTmeshes(2).X,1),1]); %Theta
% solutionomega = interpolatePGD(pgdaux,2,{[data.PGDmeshes(2).X,data.PGDmeshes(2).X]},data.OPTmeshes,data.nOfPGDterms);
% solutiontheta = interpolatePGD(pgdaux,2,{[data.PGDmeshes(2).X,data.PGDmeshes(3).X]},data.OPTmeshes,data.nOfPGDterms);

%Plot graphs
hold on
plotHandle{1} = plot(handles.axesOmega2D,2*pi./data.PGDmeshes(2).X,solution_theta);
lineHandle{1} = line([2*pi./data.snapshot{1} 2*pi./data.snapshot{1}],handles.axesOmega2D.YLim,...
    'Color',[1 0 0],'Parent',handles.axesOmega2D);
plotHandle{2} = plot(handles.axesTheta2D,180/pi.*data.PGDmeshes(3).X,solution_omega);
lineHandle{2} = line([180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2}],handles.axesTheta2D.YLim,...
    'Color',[1 0 0],'Parent',handles.axesTheta2D);
surfHandle = trisurf(delaunay(x_axis,y_axis),x_axis,y_axis,solution_omega_theta',...
    'FaceColor','interp','EdgeAlpha',0,'Parent',handles.axes3D);
rectHandle{1} = patch('XData',[2*pi./data.snapshot{1} 2*pi./data.snapshot{1} 2*pi./data.snapshot{1} 2*pi./data.snapshot{1}],...
    'YData',[handles.axes3D.YLim(1) handles.axes3D.YLim(1) handles.axes3D.YLim(2) handles.axes3D.YLim(2)],...
    'ZData',[handles.axes3D.ZLim(1) handles.axes3D.ZLim(2) handles.axes3D.ZLim(2) handles.axes3D.ZLim(1)],...
    'FaceColor','r','FaceAlpha',0.25,'EdgeAlpha',0,'Parent',handles.axes3D);
rectHandle{2} = patch('XData',[handles.axes3D.XLim(1) handles.axes3D.XLim(1) handles.axes3D.XLim(2) handles.axes3D.XLim(2)],...
    'YData',[180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2} 180/pi.*data.snapshot{2}],...
    'ZData',[handles.axes3D.ZLim(1) handles.axes3D.ZLim(2) handles.axes3D.ZLim(2) handles.axes3D.ZLim(1)],...
    'FaceColor',[1 0 0],'FaceAlpha',0.25,'EdgeAlpha',0,'Parent',handles.axes3D);
hold off

set([handles.axesOmega2D handles.axesTheta2D],'XGrid','on','YGrid','on')
set(handles.axes3D,'YLim',handles.axesTheta2D.XLim,'YDir','reverse')

%Output variable
if ~nargout
    varargout = [];
elseif nargout == 1
    varargout = {plotHandle};
elseif nargout == 2
    varargout = {plotHandle,lineHandle};
elseif nargout == 4
    varargout = {plotHandle,lineHandle,surfHandle,rectHandle};
end