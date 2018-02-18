function pointer = DrawMaxMin(axeshandle,X,variable,Name,Units)

%axeshandle = Handle to the axes in which to find max and min
%X = nodal coordinates matrix
%Variable = Variable of whose max and min we want to obtain
%Name = Variable name (to be shown in textbox)
%Units = Variable units (to be shown in textbox)

[maximum,index_max] = max(variable);
[minimum,index_min] = min(variable);

offset = 25;

%Draw plots
hold on
pointer{1}(1) = plot(axeshandle,X(index_max,1),X(index_max,2),'w.','MarkerSize',20,...
    'Visible','off','Tag','MinMaxPoint');
pointer{1}(2) = text(axeshandle,X(index_max,1)+offset,X(index_max,2),{['Max. ',Name,': ',num2str(maximum),' ',Units]},...
    'Visible','off','Tag','MinMaxPoint','BackgroundColor','w');
pointer{2}(1) = plot(axeshandle,X(index_min,1),X(index_min,2),'w.','MarkerSize',20,...
    'Visible','off','Tag','MinMaxPoint');
pointer{2}(2) = text(axeshandle,X(index_min,1)+offset,X(index_min,2),num2str(minimum),...
    'Visible','off','Tag','MinMaxPoint','BackgroundColor','w');
hold off

end