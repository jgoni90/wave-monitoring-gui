function initializeLocation(data,handles)

data = guidata(handles.MainFigure);

axes(handles.axesLocation); %Set axes as current
data.PortNames = cellstr(['Pasajes';'Mataro ']);
data.PortCoordinates = [43.323950 -1.920708;41.529646 2.444962];
data.location = imread('LOCATION.png'); % Load in a background image and display it using the correct colors
data.locationPlot = imagesc(data.location); %Load image data
Coordinates = data.PortCoordinates; %Load coordinates from data
numberofports = size(Coordinates,1); %Number of markers
normalizedCoordinates = [180 50;280 100]; %Algorithm pending
% PortXCoordinate = data.port.Coordinates(1); 
% PortYCoordinate = data.port.Coordinates(2);
hold on
for i = 1:numberofports
    PortMarkers(i) = plot(normalizedCoordinates(i,1),...
        normalizedCoordinates(i,2),'ok','Markersize',10,...
        'MarkerFaceColor','k','Visible','on','Tag',data.PortNames{i});
end
hold off

set(handles.axesLocation,'Tag','axesHandle1','visible','off');

set(PortMarkers,'buttondownfcn',{@load_pgd,handles});

data.PortMarkers = PortMarkers;
data.PlotMarkerCoordinates = normalizedCoordinates;

%Disable zoom for axesLocation and axesArrow
h = zoom;
setAllowAxesZoom(h,[handles.axesLocation handles.axesArrow handles.axesOmega2D handles.axesTheta2D],false);
hp = pan;
setAllowAxesPan(hp,[handles.axesLocation handles.axesArrow handles.axesOmega2D handles.axesTheta2D],false);
h3 = rotate3d;
setAllowAxesRotate(h3,[handles.axesArrow handles.axesLocation handles.axesOmega2D handles.axesTheta2D],false);

guidata(handles.MainFigure,data);