function Pointer = DrawMin(axes,plot,Name,Units)

x = get(plot,'Vertices',1);
y = get(plot,'Vertices',2);
z = get(plot,'Vertices',3);

[minimum,index] = min(z);

offset = 25;

hold on
Pointer(1) = plot3(axes,x(index),y(index),max(z),'r+','Markersize',20,'MarkerFaceColor',...
    'r','LineWidth',3,'Visible','on'); %Plot node point
Pointer(2) = text(x(index)+offset,y(index),max(z),...
    {[Name,': ',num2str(minimum),' ',Units],['X: ',num2str(x(index))],...
    ['Y: ',num2str(y(index))]},'Parent',axes,'BackgroundColor','w','Visible','on');
hold off

end