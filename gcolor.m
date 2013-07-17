function [h1, h2] = gcolor(varargin)
% gcolor  Plot vetical section of layered data
%
% gcolor(e)            - draw interface positions
% gcolor(h)            - draw interface positions
% gcolor(temp,e)       - shade temp and draw interface positions
% gcolor(temp,h)       - shade temp and draw interface positions
% gcolor(temp,[e|h],x) - shade temp and draw interface positions
%                        with model horizontal coordinate
%
% Shape of arguments: e(nz+1,nx), h(nz,nx), temp(nz,nx), x(nx)
%
% Grid drawing parameters:
% gcolor(...,'linear') - Linear interpolation between h-points [default]
% gcolor(...,'pcm')    - Piecewise constant (step topography)
% gcolor(...,'cplm')   - Continuous piecewise linear (non-conservative)
% gcolor(...,'plm')    - Discontinuous piecewise linear (conservative)
% gcolor(...,'quick')  - Use the orignal (fast) rendering (approximate and shifted!)
%
% Rendering of vertical representation parameters:
% gcolor(...,'vpcm'    - Use piecewise constant representation in the vertical [default]
% gcolor(...,'vplm'    - Use piecewise linear representation in the vertical
% gcolor(...,'none'    - Cheap piecewise constant only usable with shading flat
%                        (used by 'quick' option)
%
% Miscellaneous parameters:
% gcolor(...,'dots')   - Draw a black dot at the interface position
%                        This is useful for when used with 'cplm'
% gcolor(...,'nogrid') - Turns off the drawing of grid interfaces
%
% e.g. Plot a time-depth (time-series) at a single (x,y) point
% nc=netcdf('gold_output.nc');
% gcolor( nc{'temp'}(:,:,1,1)', nc{'e'}(:,:,1,1)' );
% axis([0 365 -300 0]) % zoom in on top 300 meters, first 365 days
% caxis([12 28]);colorbar % set temperature range
%
% e.g. Plot lat-depth section (single time)
% nc=netcdf('prog.nc');
% gcolor( nc{'salt'}(end,:,:,260), nc{'e'}(end,:,:,260), nc{'yh'}(:) );
% caxis([32 38]);colorbar('h') % set salt range
%
% e.g. Plot lon-depth section (single time)
% gcolor(nc{'salt'}(end,:,20,:),nc{'e'}(end,:,20,:),nc{'xh'}(:));
% caxis([33 36]);colorbar('h')
%
% Version: $Id: gcolor.m,v 1.7 2013/07/17 17:59:08 Alistair.Adcroft Exp $

rend='linear';
repr='vpcm';
draw_dots=0;
draw_grid=1;
drawMLcolors=0;
q=[];x=[];

nflds=0;
for k=1:nargin
  if isnumeric(varargin{k})
    nflds=nflds+1;
    switch nflds
      case 1,
        h=squeeze(varargin{k});
      case 2,
        q=h;
        h=squeeze(varargin{k});
      case 3,
        x=squeeze(varargin{k});
      otherwise
        error('Too many numeric arguments!')
    end % case
  elseif ischar(varargin{k})
    switch varargin{k}
      case {'quick'}
        rend='none';
        repr='none';
      case {'linear','pcm','cplm','plm','linear-orig'}
        rend=varargin{k};
      case {'none','vpcm','vplm'}
        repr=varargin{k};
      case {'dots'}
        draw_dots=1;
      case {'nogrid'}
        draw_grid=0;
      case {'noml'}
        drawMLcolors=0;
      case {'ml'}
        drawMLcolors=1;
      otherwise
        error(sprintf('Unrecognozied option: %s',varargin{k}))
    end
  else
    error(sprintf('Not sure what to do with argument %k',k))
  end % if
end % k
if isempty(q) % Only e or h was supplied so make up scalar field
  if min(h)>=0
    q=h;
  else
    q=h(2:end,:)-h(1:end-1,:);
  end
end
if isempty(x) % No coordinate supplied so make one up
  x=0.5:size(h,2);
end

ni=size(h,2);

if min(h(:))>=0
  e=[zeros(ni,1) -cumsum(h,1)']';
elseif max(h(:)) < max(abs(h(:)))
  e=h;
else
  error('Can not interpret the E|H data')
end
if size(e,1)==size(q,1) 
  repr='ints';
elseif size(e,1)~=size(q,1)+1 
  error('Size of T,[H|E] arrays are inconsistent')
end

x=x(:)';
if length(x)==size(q,2)
  xi=[1.5*x(1)-0.5*x(2) (x(2:end)+x(1:end-1))/2 1.5*x(end)-0.5*x(end-1)]; 
elseif length(x)==size(q,2)+1
  xi=x;
else
  error('X coordinate has wrong dimensions')
end

if length(xi)~=size(h,2)+1
  error('X dimension of H is inconsistent with Q or X')
end

ni=size(q,2);nk=size(q,1);nke=size(e,1);
switch lower(rend)
  case {'none'}
    X=xi; Q=q; E=e;
  case {'linear-orig'} % This is the old way that did not work with shading interp
    % Linear interpolation from cell centers
    X=zeros(1,2*ni+1); X(1:2:2*ni+1)=xi; X(2:2:2*ni)=(xi(1:ni)+xi(2:ni+1))/2;
    Q=zeros(nk,2*ni); Q(:,1:2:2*ni-1)=q; Q(:,2:2:2*ni)=q;
    E=zeros(nke,2*ni);
    E(:,3:2:2*ni-1)=(e(:,1:ni-1)+e(:,2:ni))/2;
    E(:,1)=e(:,1);
    E(:,2:2:2*ni)=e(:,:);
  case {'linear'}
    % Linear interpolation from cell centers
    X=zeros(1,3*ni+1); X(1:3:3*ni+1)=xi; X(3:3:3*ni+1)=X(4:3:3*ni+1); X(2:3:3*ni)=(xi(1:ni)+xi(2:ni+1))/2;
    Q=zeros(nk,3*ni); Q(:,1:3:3*ni-1)=q; Q(:,2:3:3*ni)=q; Q(:,3:3:3*ni)=q;
    E=zeros(nke,3*ni);
    E(:,3:3:3*ni-1)=(e(:,1:ni-1)+e(:,2:ni))/2;
    E(:,4:3:3*ni)=(e(:,1:ni-1)+e(:,2:ni))/2;
    E(:,1)=e(:,1);
    E(:,3*ni)=e(:,ni);
    E(:,2:3:3*ni)=e(:,:);
  case {'pcm'}
    % Step topography/elevations
    X=zeros(1,2*ni+1); X(1:2:2*ni+1)=xi; X(2:2:2*ni)=xi(2:ni+1);
    Q=zeros(nk,2*ni); Q(:,1:2:2*ni-1)=q; Q(:,2:2:2*ni)=q;
    E=zeros(nke,2*ni);
    E(:,1:2:2*ni-1)=e(:,:);
    E(:,2:2:2*ni)=e(:,:);
  case {'cplm'}
    % PLMish representation using mid-points for continuous edge values
    X=zeros(1,2*ni+1); X(1:2:2*ni+1)=xi; X(2:2:2*ni)=xi(2:ni+1);
    Q=zeros(nk,2*ni); Q(:,1:2:2*ni-1)=q; Q(:,2:2:2*ni)=q;
    E=zeros(nke,2*ni);
    Ee=(e(:,[1 1:ni])+e(:,[1:ni ni]))/2;
    E(:,1:2:2*ni)=Ee(:,1:ni);
    E(:,2:2:2*ni)=Ee(:,2:ni+1);
  case {'plm'}
    X=zeros(1,2*ni+1); X(1:2:2*ni+1)=xi; X(2:2:2*ni)=xi(2:ni+1);
    Q=zeros(nk,2*ni); Q(:,1:2:2*ni-1)=q; Q(:,2:2:2*ni)=q;
    E=zeros(nke,2*ni);
    sr=e(:,[2:ni ni])-e; sl=e-e(:,[1 1:ni-1]); s2=(sl+sr)/2;
    % Note: the multiplier i.f.o. sr/sl should be 2 but leads
    % to visually non-monotonic elevations. 2 is technically correct
    % for advection but 1 is required for arbitrary remapping.
    % cf. LAW's results
    s2=sign(s2).*min(min(abs(1*sr),abs(1*sl)),abs(s2));
    Er=e+0.5*s2; El=e-0.5*s2;
    E(:,1:2:2*ni)=El;
    E(:,2:2:2*ni)=Er;
  otherwise
    % Do the old ("broken") way
    X=xi;Q=q;E=e;
end
ni=size(Q,2);nk=size(Q,1);
switch lower(repr)
  case {'none'}
    % Usable with shading flat only
    EE=E;
  case {'vpcm'}
    % Use PCM representation in the vertical
    EE=zeros(2*nk+1,ni); W=Q; Q=zeros(2*nk,ni);
    EE(1:2:2*nk+1,:)=E(1:nk+1,:); EE(2:2:2*nk+1,:)=E(2:nk+1,:);
    Q(1:2:2*nk-1,:)=W;
    Q(2:2:2*nk,:)=W;
  case {'vplm'}
    % Use PLM representation in the vertical
    EE=zeros(2*nk+1,ni); W=Q; Q=zeros(2*nk,ni);
    EE(1:2:2*nk+1,:)=E(1:nk+1,:); EE(2:2:2*nk+1,:)=E(2:nk+1,:);
    St=W([2:nk nk],:)-W; Sb=W-W([1 1:nk-1],:);
    St(nk,:)=St(nk-1,:); Sb(1,:)=Sb(2,:); % Extrpolate at boundaries (NEED TO INTERIOR TOO!!!!)
    S2=(St+Sb)/2;
    S2=sign(S2).*min(min(abs(2*St),abs(2*Sb)),abs(S2));
    Q(1:2:2*nk-1,:)=W-S2/2;;
    Q(2:2:2*nk,:)=W+S2/2;
  case {'ints'}
    % Usable with shading flat only
    EE=E; EE(end+1,:)=EE(end,:);
    disp('Switching to interpolate mode from reconstruction mode. You supplied data at interfaces rather than layers.')
  otherwise
    % Use PCM representation in the vertical
    EE=E;
end

%XI=repmat(X,size(EE,1),1);
E(:,end+1)=E(:,end);
EE(:,end+1)=EE(:,end);
Q(end+1,:)=1*Q(end,:);
Q(:,end+1)=1*Q(:,end);

h1=pcolor(X,EE,Q);
if repr=='none'
error('Hi')
  shading flat
else
  shading interp
end
if draw_grid
  if drawMLcolors % assume mixed layer model
   hold on;h2=plot(X',E(6:end,:)','k','LineWidth',0.5);hold off
   hold on;h2=plot(X',E(4:5,:)','w','LineWidth',0.5);hold off
   hold on;h2=plot(X',E(1:3,:)','b','LineWidth',0.5);hold off
  else
   hold on;h2=plot(X',E(:,:)','k','LineWidth',0.5);hold off
   %hold on;h2=plot(X',E(:,:)','LineWidth',0.5);hold off
  end
end
if draw_dots
  hold on;plot((xi(1:end-1)+xi(2:end))/2,e,'k.');hold off
end
