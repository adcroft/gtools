function [T,DZ,D,zc] = gremap(eh,temp,zdz)
% Remaps layer data to z-space
% (Note: do not use this to remap interface data. Use ginterp() instead)
%
% [T,DZ]=gremap(e,temp,z)  - using interfaces, e, remap temp to z (z<=0)
% [T,DZ]=gremap(h,temp,z)  - using layer thicknesses, h, remap temp to z (z<=0)
% [T,DZ]=gremap(e,temp,dz) - using interfaces, e, remap temp to z (dz>0)
% [T,DZ]=gremap(h,temp,dz) - using layer thicknesses, h, remap temp to z (dz>0)
%
% [T,DZ]=gremap(e,temp,z) - in addition, returns DZ the actual thicknesses associated with T
% [T,DZ,D]=gremap(e,temp,z) - in addition, returns D the topography implied by the bottom interface (e)
% [T,DZ,D,zc]=gremap(e,temp,z) - in addition, returns zc, a 1D vector of notional level positions (for quick and easy plotting)
%
% Shape of arguments: e(nk+1,[ny,nx][ny][nx][1]), h(nk,[ny,nx][ny][nx][1]), temp(nk,[ny,nx][ny][nx][1]), z([NZ][NZ+1]), dz(NZ)
%   if Z(1)==0 then Z will be interpretted as interface depths to remap to (ie. to get NZ levels supply NZ+1 interfaces)
%
% e.g. 
% nc=netcdf('gold_output.nc');
% gcolor(nc{'temp'}(1,:,:,1),nc{'e'}(1,:,:,1),'pcm');
% [T,DZ]=gremap( nc{'e'}(1,:,:,:), nc{'temp'}(1,:,:,:), -[0:50:300 400:100:1000 1200:200:6000] );
% gcolor(T(:,:,1),DZ(:,:,1),'pcm');
%
% Version: $Id: gremap.m,v 1.4 2011/11/30 20:57:43 aja Exp $

tic
zdz=zdz(:);
if zdz(1)==0 % Assume interface positions
 if max(zdz)>0
   error('z(1)==0 but z(:)>0')
 end
 zi=zdz; dz=zi(1:end-1)-zi(2:end);
elseif min(zdz)>0 % Assume dz
 dz=zdz; zi=[0 -cumsum(dz)];
else % Assume mid-level depths
 zi=0;
 for k=1:length(zdz);
  dz(k)=2*(z(k)-zdz(k));
  zi(k+1)=zi(k)-dz(k);
  if k>1 & dz(k)<dz(k-1)
   disp('Warning from gremap: dz is oscillating!')
  end
 end
end
zc=zi(1:end-1)-dz/2; % For quick plotting convenience
toc

tic
tsz=size(temp);
ijsz=size(temp(:,:));
temp=temp(:,:);eh=eh(:,:);

if size(eh,1)==size(temp,1) & min(eh(:))>=0 % Assume eh is layer thicknesses
% e=[zeros(ni,1) -cumsum(eh,1)']';
  e=eh; e(end+1,:)=0; e(1,:)=0; e(2:end,:)=-cumsum(eh(:,:),1);
  h=eh;
elseif size(eh,1)==size(temp,1)+1 & max(eh(:)) < max(abs(eh(:))) % Assume eh is interface depths
  e=eh;
  h=eh(1:end-1,:)-eh(2:end,:);
else
  whos
max(abs(eh(:)))
max(eh(:))
  error('Size of temp,[H|E] arrays are inconsistent')
end
D=e(end,:);
toc

tic
reconmode='pcm';
switch lower(reconmode)
  case {'pcm'}
    % Use PCM representation in the vertical
   %Tup=temp;
   %Tdn=temp;
  case {'plm'}
    % Use PLM representation in the vertical
    error('NO PLM YET!')
  otherwise
    % Use PCM representation in the vertical
    error('Otherwise!')
end
toc

tic
% Vertical integral of temp (starting at top)
tint=0*temp;
for k=1:tsz(1)
 tint(k+1,:)=tint(k,:)+h(k,:).*(temp(k,:));
 dh=h(k,:); dh(find(dh==0))=1e-12; e(k+1,:)=e(k,:)-dh; % Should find cleaner way to do this - AJA
end
toc

tic
% Sample tint (integral of temp) on Z grid -> Tint
DZ=zeros([length(dz) ijsz(2)]);
Tint=zeros([length(zi) ijsz(2)]);
for j=1:ijsz(2)
 zstarfac=-(e(1,j)-e(end,j))./e(end,j);
 Zi=e(1,j)+zstarfac.*zi;
 Zi=max(Zi,D(j));
 DZ(:,j)=Zi(1:end-1)-Zi(2:end);
 Tint(:,j)=interp1(e(:,j),tint(:,j),Zi,'spline');
end
toc

tic
T=zeros([length(dz) ijsz(2)]);
for k=1:length(dz)
 for j=1:ijsz(2)
  Dz=DZ(k,j);
  T(k,j)=(Tint(k+1,j)-Tint(k,j))./Dz;
 end
end
T(isnan(T))=0;
toc

tic
T=reshape(T,[length(dz) tsz(2:end)]);
DZ=reshape(DZ,[length(dz) tsz(2:end)]);
D=squeeze(reshape(D,[1 tsz(2:end)]));
toc
