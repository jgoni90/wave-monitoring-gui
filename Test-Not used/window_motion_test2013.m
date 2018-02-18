function window_motion_test2013
% 
figure('WindowButtonDownFcn',@wbdcb)
ah = axes('SortMethod','childorder');
axis ([1 10 1 10])
title('Click and drag')
   function wbdcb(src,callbackdata)
     if strcmp(get(src,'SelectionType'),'normal')
        set(src,'Pointer','circle');
        cp = get(ah,'CurrentPoint');
        xinit = cp(1,1);
        yinit = cp(1,2);
        hl = line('XData',xinit,'YData',yinit,...
        'Marker','p','color','b');
        set(src,'WindowButtonMotionFcn',@wbmcb);
        set(src,'WindowButtonUpFcn',@wbucb);
     end    
 
        function wbmcb(src,callbackdata)
           cp = get(ah,'CurrentPoint');
           xdat = [xinit,cp(1,1)];
           ydat = [yinit,cp(1,2)];
           hl.XData = xdat;
           hl.YData = ydat;
           drawnow
        end
   
        function wbucb(src,callbackdata)
           if strcmp(get(src,'SelectionType'),'alt')
              set(src,'Pointer','arrow');
              set(src,'WindowButtonMotionFcn','');
              set(src,'WindowButtonUpFcn','');
           else
              return
           end
        end
  end
end