function [h] = gplananim(filename,varName,time,varargin)
% gplananim  Animate 2D plot from netcdf
%
% gplananim(filename,varname,time)
% gplananim(filename,varname,time,Z)
%
% e.g.
% >> gplananim('file.nc','salt',1:12,1);
%
% Written by A.Adcroft, Winter 2011

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

% Use global variables to control plots
global CAXIS TITLE COLORBAR PNGNAME
cax=[];
if length(CAXIS)==2
 cax=CAXIS;
end

if ischar(time)
  %error('Have not implemented character specification of time range yet!')
  t=gread(nc,'time');
  time=1:length(t);
end
if isnumeric(time)
  for t=time
   %disp(sprintf('Plotting time=%g',t))
   gplan(nc,varName,t,varargin{:});
   if isempty(cax)
     cax=caxis; % Store for future plots
   else
     caxis(cax);
   end % cax
   if length(TITLE)>0
     str=get(get(gca,'Title'),'String');
     title([TITLE ' ' str]);
   end % TITLE
   if length(COLORBAR)>0
     colorbar(COLORBAR)
   end % COLORBAR
   drawnow
   if length(PNGNAME)>0
     print('-dpng',sprintf('%s.%i5.5.png',PNGNAME,t))
   end % PNGNAME
  end % t
end

if closenc
 close(nc)
end
