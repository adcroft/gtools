function [h] = gplandiff(filename1,filename2,varName,varargin)
% gplandiff  Plot 2D variable in-x-y view from netcdf
%
% gplandiff(file1,file2,,varname,time)
% gplandiff(file1,file2,,varname,time,Z)
%
% e.g.
% >> gplandiff('file1.nc','file2.nc','salt',Inf,1);
%
% Written by A.Adcroft, Winter 2011

% Open the netcdf file
if ischar(filename1)
 closenc1=1;
 if exist(filename1,'file')
  nc1=netcdf(filename1,'nowrite');
 else
  error(['File ''' filename1 ''' does not exist'])
 end
elseif strcmp(class(filename1),'netcdf')
 closenc1=0;
 nc1=filename1;
 filename1=name(nc1);
else
 error('filename1 argument is neither a netcdf handle nor a character string')
end
if ischar(filename2)
 closenc2=1;
 if exist(filename2,'file')
  nc2=netcdf(filename2,'nowrite');
 else
  error(['File ''' filename2 ''' does not exist'])
 end
elseif strcmp(class(filename2),'netcdf')
 closenc2=0;
 nc2=filename2;
 filename2=name(nc2);
else
 error('filename1 argument is neither a netcdf handle nor a character string')
end
[q1,x,y,t,rho] = getq(nc1,varName,varargin{:});
[q2,x,y,t,rho] = getq(nc2,varName,varargin{:});

global msk
if length(msk)==length(q1)
 h=pcol(x,y,(q1-q2).*msk);
else
 h=pcol(x,y,(q1-q2));
end
xlabel('Longitude')
ylabel('Latitude')
title(['\Delta ' sprintf('%s (time=%g)',varName,t)]);

if closenc1
 close(nc1)
end
if closenc2
 close(nc2)
end

% --------------------------------------------------------------
function [q,x,y,t,rho] = getq(nc,varName,varargin)
if isempty(nc{varName}) & strcmp(varName,'vorticity')
 [u,v,xu,yu,xv,yv,t,rho] = readuv(nc,varargin{:});
 y=(yv+yv(:,[2:end end]))/2;
 x=(xu+xu([2:end end],:))/2;
 global OCEAN_GEOMETRY
 if length(OCEAN_GEOMETRY)>0
  [dxq,latq,lonq]=gread(OCEAN_GEOMETRY,'dxq',sprintf('=%g:%g',min(y),max(y)),sprintf('=%g:%g',min(x),max(x)));
  [dyq,latq,lonq]=gread(OCEAN_GEOMETRY,'dxq',sprintf('=%g:%g',min(y),max(y)),sprintf('=%g:%g',min(x),max(x)));
  dxq=dxq(1:length(yv),1:length(xu));
  dyq=dyq(1:length(yv),1:length(xu));
 else
  [yu,xv]=ndgrid(yu,xv);
  [Y,X]=ndgrid(y,x);
  Re=6370e3;
  dyq=Re*(yu([2:end end],:)-yu)*pi/180;
  dxq=Re*cos(Y*pi/180).*(xv(:,[2:end end])-xv)*pi/180;
 end
 q=(v(:,[2:end end])-v)./dxq-(u([2:end end],:)-u)./dyq; %  approximate F.D. not exact F.V. !!!
elseif isempty(nc{varName}) && ( strcmp(varName,'spd') || strcmp(varName,'ke') )
 [u,v,xu,yu,xv,yv,t,rho] = readuv(nc,varargin{:});
 y=(yv+yv([1 1:end-1],:))/2;
 x=(xu+xu(:,[1 1:end-1]))/2;
 v=(v.^2+v([1 1:end-1],:).^2)/2;
 u=(u.^2+u(:,[1 1:end-1]).^2)/2;
 if strcmp(varName,'spd')
  q=sqrt(u+v);
 else
  q=(u+v)/2;
 end
elseif isempty(nc{varName}) && ( strcmp(varName,'gspd') || strcmp(varName,'gke') )
 % Godunov style KE
 [u,v,xu,yu,xv,yv,t,rho] = readuv(nc,varargin{:});
 y=(yv+yv([1 1:end-1],:))/2;
 x=(xu+xu(:,[1 1:end-1]))/2;
 q=0*u.*u(:,[1 1:end-1]).*v.*v([1 1:end-1],:); % Allocate memory, set mask
 up=max(0,u); q=q+up.^2;
 up=min(0,u); q=q+up(:,[1 1:end-1]).^2;
 up=max(0,v); q=q+up.^2;
 up=min(0,v); q=q+up([1 1:end-1],:).^2;
 if strcmp(varName,'gspd')
  q=sqrt(q);
 else
  q=q/2;
 end
elseif isempty(nc{varName}) 
 error(['Could not find variable ''' varName ''' in file ' filename])
else
 [q,x,y,t,rho] = locgread(nc,varName,varargin{:});
end % isempty

% --------------------------------------------------------------
function [q,x,y,t,rho] = locgread(nc,varName,varargin)
rho=[];t=[];
switch length( dim(nc{varName}) )
 case {2}
  [q,y,x]=gread(nc,varName,varargin{:});
 case {3}
  [q,t,y,x]=gread(nc,varName,varargin{:});
 case {4}
  [q,t,rho,y,x]=gread(nc,varName,varargin{:});
 otherwise
  error(['Dimensions of variable,' varName ', are inconsistent with a plan view!'])
end
function [u,v,xu,yu,xv,yv,t,rho] = readuv(nc,varargin)
if ~isempty(nc{'u'})
 [u,xu,yu,t,rho] = locgread(nc,'u',varargin{:});
 [v,xv,yv,t,rho] = locgread(nc,'v',varargin{:});
elseif ~isempty(nc{'ssu'})
 [u,xu,yu,t,rho] = locgread(nc,'ssu',varargin{:});
 [v,xv,yv,t,rho] = locgread(nc,'ssv',varargin{:});
else
 error('Could not find u,v or ssu,ssv in file.')
end
