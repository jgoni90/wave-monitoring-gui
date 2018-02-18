function View_BC_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

value = get(hObject,'Value');

if strcmp(hObject.String{value},'None')
    for i = 1:size(data.boundaryid,2);
        for j = 1:length(fieldnames(data.PGDmeshes(1).elementFaceInfo));
            set(data.boundaryid{i}{j},'Visible','off')
        end
    end
elseif strcmp(hObject.String{value},'All')
    for i = 1:size(data.boundaryid,2);
        for j = 1:length(fieldnames(data.PGDmeshes(1).elementFaceInfo));
            set(data.boundaryid{i}{j},'Visible','on')
        end
    end
else
    for i = 1:size(data.boundaryid,2);
        for j = 1:length(fieldnames(data.PGDmeshes(1).elementFaceInfo));
            set(data.boundaryid{i}{j},'Visible','off')
        end
        set(data.boundaryid{i}{value-1},'Visible','on')
    end
end