function [data,varargout] = gread(filename,varname,varargin)
% gread  Read a variable from a netcdf file
%
% v = gread(filename,varname)
% v = gread(nc_handle,varname)
% [v,t,y,x] = gread(filename,varname,coord_range,coord_range,...)
% [v,t,y,x] = gread(nc_handle,varname,coord_range,coord_range,...)
%
% e.g.
% >> SST = gread('file.nc','sst');                          % Missing coordinate ranges default to ':' (i.e. all data)
% >> SST = gread('file.nc','sst',Inf);                      % Means last value in that axis
% >> SST = gread('file.nc','sst',:,:);                      % : means all data on that axis
% >> SST = gread('file.nc','sst',:,5:10,:);                 % 5:10 means indices 5 to 10
% >> SST = gread('file.nc','sst',:,'yh=-5:5','xh=10:35');   % Coordinate names must match the netcdf dimensions
% >> SST = gread('file.nc','sst',:,'-5:5','10:35');         % Assumes order of arguments match netcdf dimensions
% >> SST = gread(nc,'sst');                                 % nc is a netcdf handle obtained with nc=netcdf(filename)
%
% Written by A.Adcroft, Fall 2011

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

% Default ranges
for j=1:length(variable_dimensions);
 drange{j}=':';
end

% For each argument, adjust the relevent coordinate range
nargs=nargin-2;
for j=1:nargs
 drange=readargs(varargin{j},j,drange,nc,variable_dimensions);
end

% Read data
try
 cmd=['data=nc{varname}(' add_commas(drange) ');'];
 cmd=['data=v_nc(' add_commas(drange) ');'];
 eval(cmd);
 if verbose~=0
   disp(['Reading with nc{''' varname '''}(' add_commas(drange) ');'])
 end
catch ME
 disp('An error occured reading the data. The netcdf command attempted was:')
 disp(['  ' cmd])
 disp(['where varname=' varname ' which has dimensions:'])
 disp(ncnames(variable_dimensions))
 error('gread failed to read any netcdf data.')
end

% Turn missing data into NaN's
fv=fillval(v_nc);
if isempty(fv)
 fv=v_nc.missing_value(:);
end
if isempty(fv)
 fv=nc.missing_value(:);
end
if ~isempty(fv)
 data(data==fv)=NaN;
end
 
% For each coordinate output argument, read the coordinate data
for j=1:(nargout-1)
 eval(['varargout{j}=nc{name(variable_dimensions{j})}(' drange{j} ');'])
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
function [drange] = readargs(arg,argpos,drange,nc,variable_dimensions)
if ischar(arg)
 if isempty(regexp(arg,'=')) % No variable name, assume regular netcdf syntax
% drange{argpos}=arg
  if arg==':'
    drange{argpos}=':';
  else
   coord=nc{name(variable_dimensions{argpos})}(:);
   drange{argpos}=coords2string(coord,arg);
  end % :
 else
  toks=regexp(arg,'=','split');
  if length(toks)==2 && length(toks{1})==0 && length(toks{2})>0 % '=range' syntax
   toks{1}=name(variable_dimensions{argpos});
  elseif length(toks)~=2 || length(toks{1})==0 || length(toks{2})==0
   error(['Unknown syntax ''' arg ''''])
  end
  newargpos=get_dim_pos(variable_dimensions,toks{1});
  coord=nc{toks{1}}(:);
  drange{newargpos}=coords2string(coord,toks{2});
 end % =
else % arg is a number, convert it to a string
 if length(arg)==1
  if isinf(arg)
   drange{argpos}='end';
  else
   drange{argpos}=sprintf('%i',arg);
  end % Inf
 else
  drange{argpos}=['[' sprintf('%g ',arg) ']'];
 end % length
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [str] = coords2string(coord,expr)
[is,ie]=coords2indices(coord,expr);
if isempty(ie) || is==ie
 str=sprintf('%i',is);
else
 str=sprintf('%i:%i',is,ie);
end % isempty(ie)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [is,ie] = coords2indices(coord,expr)
rng=sscanf(expr,'%f:%f');
is=find(coord(2:end)>rng(1) & coord(1:end-1)<=rng(1));
if isempty(is)
 is=1;
end
if length(rng)==1
 ie=is;
else
 ie=find(coord(2:end)>rng(2) & coord(1:end-1)<=rng(2))+1;
end % length(rng)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [dpos] = get_dim_pos(dimensions,dimname)
dpos=[];
for j=1:length(dimensions)
 if strcmp(name(dimensions{j}),dimname)
  dpos=j;
  break
 end
end
if isempty(dpos)
 disp('Valid dimensions for variable are:')
 disp(ncnames(dimensions))
 error(['Could not find dimension named ''' dimname '''']);
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [str] = add_commas(cells)
str=cells{1};
for j=2:length(cells)
 str=[str ',' cells{j}];
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
