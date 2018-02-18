function UpdatePoint(plot,pointer,NewTab,Name,Units)

Variable = get(plot{NewTab}(1),'FaceVertexCData');
set(pointer(2),'String',{[Name,': ',num2str(Variable(pointer(3))),' ',Units],...
    ['X: ',num2str(get(pointer(1),'XData'))],...
    ['Y: ',num2str(get(pointer(1),'YData'))]});

uistack(pointer,'top')

end