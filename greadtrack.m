function [data,varargout] = greadtrack(filename,varname,track,varargin)
% gread  Read a variable from a netcdf file
%
% v = greadtrack(filename,varname,track)
% v = greadtrack(nc_handle,varname,track)
% [v,t,x] = gread(filename,varname,track,coord_range,coord_range,...)
% [v,t,z,x] = gread(nc_handle,varname,track,coord_range,coord_range,...)
%
% e.g.
% >> track.x=[-210 -100 -100];
% >> track.y=[-60 -60 10];
% >> SST = gread('file.nc','sst',track);           % Missing coordinate ranges default to ':' (i.e. all data)
% >> SST = gread('file.nc','sst',track,Inf);       % Means last value in that axis
% >> SST = gread('file.nc','sst',track,:,:);       % : means all data on that axis
% >> SST = gread('file.nc','sst',track,:,5:10);    % 5:10 means indices 5 to 10
% >> SST = gread(nc,'sst',track);                  % nc is a netcdf handle obtained with nc=netcdf(filename)
%
% Written by A.Adcroft, Spring 2012

global verbose

% Open the netcdf file
if ischar(filename)
 closenc=1;
 if exist(filename,'file')
  nc=netcdf(filename,'nowrite');
 else
  error(['File ''' filename ''' does not exist'])
 end
elseif strcmp(class(filename),'netcdf')
 closenc=0;
 nc=filename;
 filename=name(nc);
else
 error('filename argument is neither a netcdf handle nor a character string')
end

try
 v_nc=get_var_handle(nc,varname);
catch
 disp(['Error: variable ''' varname ''' was not found. Possible variables are:'])
 disp(ncnames(var(nc)))
 error('Could not find variable in netcdf file')
end
global_attributes=att(nc);
global_dimensions=dim(nc);
variable_attributes=att(v_nc);
variable_dimensions=dim(v_nc);

if isfield(track,'x') & isfield(track,'y')
 if length(track.x) ~= length(track.y)
  error('Track coord vectors must be the same length')
 end
 if length(track.x)<2
  error('Track coord vectors must have at least two points')
 end
   
 miny=min(track.y); maxy=max(track.y);
 yrng=sprintf('=%g:%g',miny,maxy);
 minx=min(track.x); maxx=max(track.x);
 xrng=sprintf('=%g:%g',minx,maxx);
else
 error('Have not implemented i,j tracks yet')
end

switch length(varargin)
 case {0} % No space/time
  [data_block,y,x]=gread(nc,varname,yrng,xrng);
  varargout={};
 case {1}
  [data_block,t,y,x]=gread(nc,varname,varargin{:},yrng,xrng);
  varargout{1}=t;
 case {2}
  [data_block,t,z,y,x]=gread(nc,varname,varargin{:},yrng,xrng);
  varargout{1}=t;
  varargout{2}=z;
 otherwise
  error('Not coded for 4+ varargin{:}')
end

% Track coordinates
[xf,yf,sf] = generate_track(x,y,track.x,track.y);
varargout{end+1}=sf;
varargout{end+1}=yf;
varargout{end+1}=xf;

% Interpolate data to track
switch ndims(data_block)
 case {2} % No space/time
  data=interp2(x,y,data_block,xf,yf,'nearest');
 case {3} % One space/time
  data=zeros(size(data_block,1),length(sf));
  for k=1:size(data_block,1)
   data(k,:)=interp2(x,y,squeeze(data_block(k,:,:)),xf,yf,'nearest');
  end % k
 case {4} % Two space/time
  data=zeros(size(data_block,1),size(data_block,2),length(sf));
  for n=1:size(data_block,1)
   for k=1:size(data_block,2)
    data(n,k,:)=interp2(x,y,squeeze(data_block(n,k,:,:)),xf,yf,'nearest');
   end % k
  end % n
 otherwise
  error('Not coded for this size data_block')
end

if closenc
 close(nc)
end

% ----------------------------------------------------------------------------
function [vh] = get_var_handle(nc,varname)
vh=nc{varname};
if ~isempty(vh)
 return
end
% Try a different case
for v = var(nc)
 if strcmp( lower(varname), lower(name(v{1})) )
  vh=nc{name(v{1})};
 end
end
if isempty(vh)
 error(['Could not find variable named ''' varname '''']);
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [xf,yf,s] = generate_track(xd,yd,xt,yt)
res=3;
% Resolution of data grid
if length(xd)>1
 dx=min(diff(xd));
else
 dx=1;
end
if length(yd)>1
 dy=min(diff(yd));
else
 dy=1;
end
% Interpolate coordinates to track
npnts=length(xt);
is=1;
for k=1:npnts-1
 sdx=xt(k+1)-xt(k);
 sdy=yt(k+1)-yt(k);
%ni=ceil(max(sdx/dx,sdy/dy)); ni=max(ni,2);
 ni=sqrt((sdx/dx)^2+(sdy/dy)^2); ni=res*max(ni,2);
 is=[is k+(1/ni:1/ni:1)];
end
xf=interp1(1:npnts,xt,is);
yf=interp1(1:npnts,yt,is);
% Along track distance
Re=6370e3;
dx=diff(xf)*pi/180*Re.*cos(pi/180*(yf(1:end-1)+yf(2:end))/2);
dy=diff(yf)*pi/180*Re;
ds=sqrt(dx.^2+dy.^2);
s=cumsum([0 ds]);
s=s-s(1); % km from beginning of track
s=s/Re*180/pi; % Convert to degrees latitude
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
