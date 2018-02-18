function Pointer = Draw3DPoint(axes,xyz,plot,Name,Units)

%axes = Axes in which to plot Pointer
%xyz = CurrentPoint xyz coordinates
%Name = Variable name, to be displayed in textbox
%Units = Variable units, to be displayed in textbox

offset = 25;

hold on
Pointer(1) = plot3(axes,xyz(1),xyz(2),max(plot.Vertices(:,3)),'w+',...
    'Markersize',20,'MarkerFaceColor','w','LineWidth',3,'Visible','on'); %Plot node point
Pointer(2) = text(xyz(1)+offset,xyz(2),max(plot.Vertices(:,3)),...
    {[Name,': ',num2str(xyz(3)),' ',Units],['X: ',num2str(xyz(1))],...
    ['Y: ',num2str(xyz(2))]},'FontUnits','pixels','FontSize',16,'Parent',axes,'BackgroundColor','w','Visible','on');
hold off

end