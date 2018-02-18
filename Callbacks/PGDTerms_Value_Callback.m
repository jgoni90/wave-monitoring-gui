function PGDTerms_Value_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

%New # pgd terms
value = str2double(get(hObject,'String'));
maxterms = size(data.pgd.RB{1},2);
if value > maxterms
    value = maxterms;
end
set(hObject,'String',num2str(value))
data.nOfPGDterms = value;

%Draw solution
DoPlot(data,handles,true);

guidata(handles.MainFigure,data);