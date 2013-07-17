function [S] = gstruct(filename,readall)
% gstruct  Read a netcdf file into a structure
%
% [S] = gstruct(filename,readall)
%
% e.g.
% >> S = gstruct('file.nc');
% >> S = gstruct('file.nc',1);  % Reads all data
% >> S = gstruct(nc);           % Uses handle from nc=netcdf(filename)
%
% Written by A.Adcroft, Fall 2011

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

% Creat an empty structure
S.filename=filename;
S.attributes=[];
S.dimensions=[];

if ~exist('readall','var')
  readall=0;
end

% For each global attribute in the netcdf file
S.attributes=getatts(att(nc));

% For each dimension in the netcdf file
d=dim(nc);
dimensions=[];
recdims={};
for n=1:length(d)
  dimensions=setfield(dimensions,name(d{n}),d{n}(:));
  if isrecdim(d{n})
   recdims{end+1}=name(d{n});
  end
end % n
dimensions.record_dimensions=recdims;
S.dimensions=dimensions;

% For each variable in the netcdf file
v=var(nc);
data=[];
for n=1:length(v)
  V=[];
  % For each attribute of this variable
  attributes=getatts(att(v{n}));
  if ~isempty(attributes)
    V.attributes=attributes;
  end
  V.dimensions=ncnames(dim(v{n}));
  V.datatype=datatype(v{n});
  S=setfield(S,name(v{n}),V);
  if readall
   data=setfield(data,name(v{n}),v{n}(:));
  end
end % n
if closenc
 close(nc)
end
if readall
 S.data=data;
end

% ----------------------------------------------------------------------------
function [A] = getatts(a)
A=[];
for n=1:length(a)
 A=setfield(A,name(a{n}),a{n}(:));
end
