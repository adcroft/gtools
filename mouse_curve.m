function [track] = mouse_curve()

disp('Left click to add point. Middle click to delete last point. Right click to end')
but=0;
x=[];y=[];
hold on;h=plot([NaN x],[NaN y],'k.-');hold off
while(but<3)
  [xb,yb,but] = ginput(1);
  switch but
    case 1
      x(end+1)=xb;
      y(end+1)=yb;
    case 2
      if length(x)>0
        x(end)=[]; y(end)=[];
      end
    case 3
      break
    otherwise
      error('Unknown mouse button!')
  end
  set(h,'XData',x)
  set(h,'YData',y)
end
delete(h)

track.x=x;
track.y=y;
