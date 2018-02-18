function plotid = DoArrow(data,handles)

%Arrow definition points
pA = [pi -0.6 ; 0 0];

%Wave function
wp = data.snapshot{1}*data.arrowdata.a1 + data.arrowdata.b1;
x = 0:0.01:pi;
y = sin(wp*x);

%Rotation matrix
theta = data.snapshot{2} - pi;
rotmat = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];

%Rotated points
xyr = rotmat * [x ; y];
pAr = rotmat * pA;

%Plot arrow, wave function and selection point
axes(handles.axesArrow)
plotid(1) = arrow(pAr(:,1),pAr(:,2),'Width',3,'BaseAngle',45,'tipAngle',30);
hold on
plotid(2) = plot(xyr(1,:),xyr(2,:),'linewidth',2.5,'color','b');
hold on
if data.userSelection, c = 'r'; else c = 'k'; end
plotid(3) = plot(data.arrowdata.p0(1),data.arrowdata.p0(2),...
    [c 'o'],'markerfacecolor',c,'markersize',13);
uistack(plotid(3),'top')