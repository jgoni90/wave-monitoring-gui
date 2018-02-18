function [] = gui_which_axes()
% Display which axes the pointer is over....
S.fh = figure('units','pixels',...
              'position',[560 228 560 420],...
              'menubar','none',...
              'name','gui_which_axes',...
              'numbertitle','off',...
              'resize','off');
% Now make a simple plot.
x = 0:.1:2*pi;
S.ax(1) = subplot(2,2,1);
plot(x,sin(x))
S.ax(2) = subplot(2,2,2);
plot(x,cos(x))
S.ax(3) = subplot(2,2,3);
plot(x,tan(x))
S.ax(4) = subplot(2,2,4);
plot(x,x.^2)
set(S.ax,'unit','pix');
% Fill the structure with data.
S.AXP = get(S.ax,'pos');
S.tx(1) = uicontrol('style','tex',...
                    'unit','pix',...
                    'posit',[50 395 250 22],...
                    'backg',get(S.fh,'color'),...
                    'fontsize',12,'fontweight','bold',... 
                    'string','Current Pointer Axes:');
% This textbox will display the current position of the mouse.
S.tx(2) = uicontrol('style','tex',...
                    'unit','pix',...
                    'position',[310 395 120 25],...
                    'backg',get(S.fh,'color'),...
                    'fontsize',12,'fontweight','bold' );
set(S.fh,'windowbuttonmotionfcn',{@fh_wbmfcn,S}) % Set the motion detector.
function [] = fh_wbmfcn(varargin)
% WindowButtonMotionFcn for the figure.
S = varargin{3};  % Get the structure.
set(S.tx(2),'string','None')
F = get(S.fh,'currentpoint');  % The current point w.r.t the figure.
% Figure out of the current point is over the axes or not -> logicals.
for ii = 1:4
    tf1 = S.AXP{ii}(1) <= F(1) && F(1) <= S.AXP{ii}(1) + S.AXP{ii}(3);
    tf2 = S.AXP{ii}(2) <= F(2) && F(2) <= S.AXP{ii}(2) + S.AXP{ii}(4);
      if tf1 && tf2
          set(S.tx(2),'str',['Axes ' num2str(ii)])
          break
      end
  end