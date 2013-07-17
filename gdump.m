function [] = gstruct(filename)
% gdump  Concise alternative to ncdump
%
% gstruct(filename)
% gstruct(nc)
%
% e.g.
% >> gdump('file.nc');
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
 
% For each global attribute in the netcdf file
%S.attributes=getatts(att(nc));

disp('Dimensions')
disp('==========')
global_dims=dim(nc);
str='';
for d=1:length(global_dims)
 if isrecdim(global_dims{d})
  str=[str sprintf('*%s* ',name(global_dims{d}))];
 else
  str=[str sprintf('%s ',name(global_dims{d}))];
 end
end
disp(str)

% For each coordinate variable in the netcdf file
disp('Coordinates')
disp('===========')
v=var(nc);
for n=1:length(v)
  isCoordinate=0;
  for g=1:length(global_dims)
   if strcmp(name(global_dims{g}),name(v{n}))
    isCoordinate=1;
   end
  end % d
  if isCoordinate==0
   continue
  end
  dims=dim(v{n});
  if length(v{n})>2
   disp([sprintf('%s(%i)',name(v{n}),length(v{n})) sprintf('=[%g, %g, ... %g]',v{n}(1),v{n}(2),v{n}(end))])
  else
   disp([sprintf('%s(%i)',name(v{n}),length(v{n})) sprintf('=[%g ... %g]',v{n}(1),v{n}(end))])
  end
end % n
% For each variable in the netcdf file
disp('Variables')
disp('=========')
v=var(nc);
for n=1:length(v)
  isVariable=1;
  for g=1:length(global_dims)
   if strcmp(name(global_dims{g}),name(v{n}))
    isVariable=0;
   end
  end % d
  if isVariable==0
   continue
  end
  dims=dim(v{n});
  fprintf('%s(',name(v{n}))
  for d=1:length(dims)
   fprintf('%s[%i]',name(dims{d}),length(dims{d}))
   if d<length(dims)
    fprintf(',')
   end
  end % d
  fprintf(')\n')
end % n
if closenc
 close(nc)
end

% ----------------------------------------------------------------------------
function [A] = getatts(a)
A=[];
for n=1:length(a)
 A=setfield(A,name(a{n}),a{n}(:));
end
