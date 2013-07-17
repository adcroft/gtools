function [list] = grank3(nCores,Grid,nPEs)
% grank3  Lists all possible decompositions upto PE count
%
% grank3(nCores,Grid,nPEs)
%
% Lists all possible decomposition using from 1 to nPEs for a given
% model grid size (Grid) and cores per node (nCores).
%
%  nCores - (scalar) number of cores per node/blade
%  Grid   - (2-vector) number of points in model grid
%  nPEs   - number of PEs to use
%
% See also grank2.
%
% e.g. 24-core, 360x210 grid, upto 96 PEs
% grank3(24,[360 210],96)
%
% To find unit-aspect ratio points/core
% q=grank3(24,[1440 1080],2400);
% a=q(:,6)./q(:,7); j=find( abs(a-1)<.1 )
% q(j,:)
%
% To find unit-aspect ratio points/core AND core-layout is 4x6 or 6x4
% q=grank3(24,[1440 1080],2400);
% a=q(:,6)./q(:,7); b=q(:,4)-q(:,5); j=find( abs(a-1)<.1 & abs(b)<3 )
% q(j,:)


if nargout~=0
 list=[];
end
npe=nPEs;
for b=1:nPEs/nCores
 if nargout==0
  grank2(nCores,Grid,b*nCores);
 else
  this=grank2(nCores,Grid,b*nCores);
 end
 if exist('this','var')
  for j=1:size(this,1)
   list(end+1,:)=this(j,:);
  end
 end
end
