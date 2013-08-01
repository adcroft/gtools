function [rgb] = interpcolormap(varargin)
% interpcolormap  Creates fancy colormaps simply
%
% interpcolormap(rgb);
% interpcolormap(x,rgb);
% interpcolormap(dx,rgb);
% interpcolormap(dx,rgb,stretching); % 'l','h','ll','hh','c'
% interpcolormap(name); % 'bwr','jet','gbr'
% interpcolormap(name,n); % 'bwr','jet','gbr'
%
% interpcolormap([[0 .2 0];[0 1 0];[1 1 0];[0 1 1];[0 0 1];[1 0 1];[1 0 0];[1 .8 .8]]);
% interpcolormap([[0 0 .1];[0 0 1];[0 1 1];[0 1 0];[1 1 0];[1 0 1];[1 0 0];[1 .9 .8]],[1 1 1 1 1.3 .7 1.3]);
% interpcolormap([[0 0 .1];[0 0 1];[0 1 1];[0 1 0];[1 1 0];[1 0 1];[1 0 0];[1 .9 .8]],[1 1 1 1 1.3 .7 1.3],'l');
% interpcolormap([[0 0 .1];[0 0 1];[0 1 1];[0 1 0];[1 1 1];[1 1 0];[1 0 1];[1 0 0];[.3 0 0]],[1 1 1 .2 .2 1.3 .7 1]); % Blue-white-red
% cm=interpcolormap([[0 0 0];[1 1 1]]);
% interpcolormap('bwr')
% interpcolormap('gbr')
%
% interpcolormap('show')	% Display colormap graphically
% interpcolormap('show','jet')	% Display colormap graphically
%
% Available maps:
%  'bwr' 'rwb' 'wbwrw' 'bgr' 'lots' 'gbr' 'bgypr' 'jet' 'land' 'water' 'water_beach'
% Available modifiers:
%  'l' 'h' 'll' 'hh' 'c' 'cc'
%
% Written by A. Adcroft, Princeton University, 2011.


if nargin==0
  error('Arguments?')
end
% Defaults
n=size(colormap,1);
RGB=[[0 .2 0];[0 1 0];[1 1 0];[0 1 1];[0 0 1];[1 0 1];[1 0 0];[1 .8 .8]];
x=x_from_dx([1 1 1 1 1.3 .7 1.3]);
pow=1;
show=0;
reverse=0;
x0=0;
for j=1:nargin
 if ischar(varargin{j})
  switch lower(varargin{j})
   case 'show'
    show=1-show;
   case 'reverse'
    reverse=1-reverse;
   case 'l'
    pow=2; x0=0;
   case 'h'
    pow=2; x0=1;
   case 'll'
    pow=3; x0=0;
   case 'hh'
    pow=3; x0=1;
   case 'c'
    pow=1.2; x0=0.5;
   case 'cc'
    pow=1.5; x0=0.5;
   case 'w'
    pow=0.7; x0=0.5;
   case 'ww'
    pow=0.5; x0=0.5;
   case 'bwr'
    RGB=[[0 0 .1];[0 .2 1];[1 1 1];[1 .2 0];[.3 0 0]];
    x=x_from_dx([1 .4 .4 1]);
   case 'rwb'
    RGB=[[0.3 0 0];[1 .2 0];[1 1 1];[0 .2 1];[0 0 .1]];
    x=x_from_dx([1 .4 .4 1]);
   case 'bgwpr'
    RGB=[[0 0 .1];[0 0 1];[0 1 1];[0 1 0];[1 1 1];[1 1 0];[1 0 1];[1 0 0];[.3 0 0]];
    x=x_from_dx([1 1 1 .2 .2 1.3 .7 1]);
   case 'wbwrw'
    RGB=[[.9 .9 1];[0 0 .2];[0 0 1];[0 1 1];[0 1 0];[1 1 1];[1 1 0];[1 0 1];[1 0 0];[.3 0 0];[1 .9 .9]];
    x=x_from_dx([.75 .25 1 1 .2 .2 1.3 .7 .25 .75]);
   case 'bgr'
    RGB=[[0 0 .1];[0 0 1];[0 1 1];[0 .5 0];[0 1 0];[1 1 0];[1 0 1];[.4 0 .4];[1 0 0];[1 0.9 0.9]];
    x=x_from_dx([1 1 1 1 1 1 1 1 1]);
   case 'bgypr'
    RGB=[[0 0 .1];[0 0 1];[0 1 1];[0 .5 0];[0 1 0];[1 1 0];[1 0 1];[.4 0 .4];[1 0 0];[1 0.9 0.9]];
    x=x_from_dx([1 .8 1 .5 .2 .5 1 1 1]);
   case 'lots'
    RGB=[[0 0 .1];[0 0 1];[0 1 1];[0 .5 0];[0 1 0];[.5 .5 0];[1 1 0];[1 0 1];[.4 0 .4];[1 0 0];[.5 0 0];[1 0.9 0.7]];
    x=x_from_dx([.7 1 1 1 1 1 1 1 1 1 1]);
   case 'gbr'
    RGB=[[0 .2 0];[0 1 0];[1 1 0];[0 1 1];[0 0 1];[1 0 1];[1 0 0];[1 .8 .8]];
    x=x_from_dx([1 1 1 1 1.3 .7 1.3]);
   case 'jet'
    RGB=[[0 0 .5];[0 0 1];[0 1 1];[1 1 0];[1 0 0];[.5 0 0]];
    x=x_from_dx([0.5 1 1 1 0.5]);
   case 'rjet'
    RGB=[[.5 0 0];[1 0 0];[1 1 0];[0 1 1];[0 0 1];[0 0 .5]];
    x=x_from_dx([0.5 1 1 1 0.5]);
   case 'land'
    RGB=[[0 0.2 0];[0 .7 0];[0 1 .3];[1 1 .3];[.7 .5 .5]];
    x=x_from_dx([.7 .3 .5 1]);
   case 'water'
    RGB=[[0 0 .1];[.3 .3 .4];[.3 .5 .6];[.1 .2 .3];[.3 .3 .7];[.8 .8 1]];
    x=x_from_dx([1 1 1 1 1]);
   case 'water_beach'
    RGB=[[0 0 .1];[.3 .3 .4];[.3 .5 .6];[.1 .2 .3];[.3 .3 .7];[.8 .8 1];[.9 .9 0]];
    x=x_from_dx([1 1 1 1 1 1e-6]);
   case 'land_water'
    if nargin==2
      Hrng=varargin{2};
      nw=round(-n*Hrng(1)/diff(Hrng))
      nl=n-nw
      cmw=interpcolormap('water',nw);
      cml=interpcolormap('land',nl);
      rgb=[cmw;cml];
      colormap(rgb)
      whos
      return
    else
      error('Oops')
    end
   otherwise
    error('Unrecognized string')
  end % switch
 elseif numel(varargin{j})==1
  n=varargin{j};
 elseif size(varargin{j},2)==3 & size(varargin{j},1)>1 % colormap
  RGB=varargin{j};
  m=size(RGB,1);
  x=(0:m-1)/(m-1);
 elseif prod(size(varargin{j}))==numel(varargin{j}) % vector
  dx=varargin{j};
  if length(dx)==size(RGB,1)
   x=cumsum([0 dx]); x=(x(1:end-1)+x(2:end))/2; x=(x-min(x))/(max(x)-min(x));
  elseif length(dx)==size(RGB,1)-1
   x=cumsum([0 dx]); x=(x-min(x))/(max(x)-min(x));
  end
 end
end % n

xc=x-x0; x=abs(xc).^pow.*sign(xc)+x0; x(1)=0; x(end)=1;
X=0:1/(n-1):1;
rgb=interp1(x,RGB,X);
if reverse
 rgb=rgb(end:-1:1,:);
end
colormap(rgb);
if nargout==0
 clear rgb
end
if show
 showcolormap;
end


function [x] = x_from_dx(dx)
x=cumsum([0 dx]); x=(x-min(x))/(max(x)-min(x));

function [] = showcolormap
subplot(211)
cm=colormap;
pcol(1:size(cm,1));

x=1:size(cm,1);
subplot(212)
%plot(cm(:,[3 2 1]))
plot(x,cm(:,3),'bx-',x,cm(:,2),'go-',x,cm(:,1),'r.-')
axis([1 size(cm,1) 0 1])
legend('b','g','r')
