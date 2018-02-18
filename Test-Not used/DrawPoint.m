function Pointer = DrawPoint(axes,xy,X,Variable,Name,Units)

%Vertex = CurrentPoint
%X = Vertex spatial distribution matrix
%Variable = Variable of current plot

IDX = knnsearch(X,xy);
Up = Variable(IDX);
cpc = [X(IDX,1) X(IDX,2)];
offset = 25;

hold on
Pointer(1) = plot3(axes,cpc(1),cpc(2),Up+1,'w+','Markersize',20,'MarkerFaceColor',...
    'w','LineWidth',3,'Visible','on'); %Plot node point
Pointer(2) = text(cpc(1)+offset,cpc(2),Up+1,...
    {[Name,': ',num2str(Up),' ',Units],['X: ',num2str(cpc(1))],...
    ['Y: ',num2str(cpc(2))]},'Parent',axes,'BackgroundColor','w','Visible','on');
hold off

end