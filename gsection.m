function [h] = gsection(filename,varName,Time,Y,X,efilename)
% gsection  Plot section of layered variable from netcdf
%
% gsection(filename,varname,time,Y,X)
% gsection(filename,varname,time,Y,X,filename_for_e_or_h)
%
% e.g.
% >> gsection('file.nc','salt',Inf,'=-30',:);
%
% Written by A.Adcroft, Fall 2011

% Open the netcdf file
[nc,filename,closenc] = opennetcdf(filename);

if exist('efilename','var')
 [enc,efilename,eclosenc] = opennetcdf(efilename);
else
 enc=nc;
 eclosenc=0;
end
if ~isempty( enc{'e'} )
 enm='e';
elseif ~isempty( enc{'E'} )
 enm='E';
elseif ~isempty( enc{'h'} )
 enm='h';
elseif ~isempty( enc{'H'} )
 enm='H';
else
 error('Can not find "e" or "h" as variables in the netcdf file')
end
[e,t,rho,y,x]=gread(enc,enm,Time,':',Y,X);
if length(y)==1 & length(x)==1
 error('One of either X or Y must be multi-valued')
elseif length(y)>1 & length(x)>1
 error('Either X or Y must be single valued')
end
if length(t)>1
 error('Time must be single valued')
end
[q,t,rho,y,x]=gread(nc,varName,Time,':',Y,X);

% Ignore interface hieghts when bottom is >=DMIN
global DMIN
if ~isempty(DMIN)
  if min(e(:))<0
   j=find(e(end,:)>=DMIN)
   q(:,j)=NaN;
   e(:,j)=0;
  end
end

if length(x)==1
 h=gcolor(q,e,y);
 xlabel('Y');
else
 h=gcolor(q,e,x);
 xlabel('X');
end
ylabel('Depth (m)')
title(sprintf('%s (time=%g)',  regexprep(varName,'_','\\_'), t));

if closenc
 close(nc)
end
if eclosenc
 close(enc)
end

%---------------------------------------------------------------
function [nc,filename,closenc] = opennetcdf(filename);
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
