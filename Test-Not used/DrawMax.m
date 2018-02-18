function Pointer = DrawMax(axes,plot,Name,Units)

%axes = Axes in which to plot Pointer
%Name = Variable name, to be displayed in textbox
%Units = Variable units, to be displayed in textbox

x = plot.Vertices(:,1);
y = plot.Vertices(:,2);
z = plot.Vertices(:,3);

[maximum,index] = max(z);

offset = 25;

hold on
Pointer(1) = plot3(axes,x(index),y(index),maximum,'r+','Markersize',20,'MarkerFaceColor',...
    'r','LineWidth',3,'Visible','on'); %Plot node point
Pointer(2) = text(x(index)+offset,y(index),maximum,...
    {[Name,': ',num2str(maximum),' ',Units],['X: ',num2str(x(index))],...
    ['Y: ',num2str(y(index))]},'Parent',axes,'BackgroundColor','w','Visible','on');
hold off

end