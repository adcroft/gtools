function [varargout] = stats(A,varargin)
% stats  Reporst the basic statistics of an array
%
% stats(A) writes the basic statistics of A to the terminal which
% are i) the minimum finite value
%    ii) the maximum finite value
%   iii) the mean of the finite values
%    iv) the S.D. of the finite values ( RMS of [A-mean] )
%     v) the fraction of non-finite elements excluded from calculations
%
% e.g.
%  >> stats(topo)
% Min -4555     Max 0         Mean -2331.07  SD 1207.3441   N-Z 1024/1024
%
% [Min Max Mean SD]=stats(topo); returns the statistics in Min, Max,
% Mean and SD and does not write to the terminal.

% $Header: /home/aja/CVS/gtools/stats.m,v 1.2 2011/12/02 20:06:21 aja Exp $
% $Name:  $

if (isempty(A))
 A=NaN;
end

A=A(:);
switch (nargin)
 case {1}
  label='';
 case {2}
  label=varargin{1};
 otherwise
  error('Too many input arguments given')
end
 
%ii=find(A~=0);
ii=find(isfinite(A));
if isempty(ii)
 ii=1;
end
sZ=prod(size(A));
minA=min(A(ii));
maxA=max(A(ii));
meanA=mean(A(ii));
sdA=sqrt( mean( (A(ii)-meanA).^2 ) );
rmsA=sqrt( mean( A(ii).^2 ) );
nZ=sum(~isfinite(A));
switch max(nargout)
 case {0}
 case {1}
  varargout(1)={minA};
 case {2}
  varargout(1)={minA};
  varargout(2)={maxA};
 case {3}
  varargout(1)={minA};
  varargout(2)={maxA};
  varargout(3)={meanA};
 case {4}
  varargout(1)={minA};
  varargout(2)={maxA};
  varargout(3)={meanA};
  varargout(4)={sdA};
 otherwise
  error('Too many return arguments requested')
end

if (nargout==0)
% disp( ...
%  sprintf('Min %0.5g  Max %0.5g  Mean %0.5g  SD %0.5g   NaN %i/%i %s',...
%        minA,maxA,meanA,sdA,nZ,sZ,label) );
  disp( ...
   sprintf('Min %0.5g  Max %0.5g  Mean %0.5g  RMS %0.5g   NaN %i/%i %s',...
         minA,maxA,meanA,rmsA,nZ,sZ,label) );
end
