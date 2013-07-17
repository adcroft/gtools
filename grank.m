function [] = grank(Sc,Sd,Sn,offset,lims)
% grank   Create MPICH_RANK_REORDER for given composition
%
% grank(Sc,Sd,Sn)
% grank(Sc,Sd,Sn,offset)
% grank(Sc,Sd,Sn,offset,lims)
%
% Create a MPICH_RANK_REORDER file for use on a heirarchically
% connect HPC
%
%  Sc is a 2-vector of the layout of cores on a chip
%  Sd is a 2-vector of the layout of chips on a node
%  Sn is a 2-vector of the layout of nodes
%
% Optional:
%  offset - is the starting MPI rank for the first core
%
% e.g.
% 24 cores on 1 node
% >> grank([2 3],[2 2],[1 1])  
% 24 cores on 1 node, starting at 72
% >> grank([2 3],[2 2],[1 1],72)  
% 48 cores on 1 node, using 2 threads (effectively halves the cores per node)
% >> grank([1 3],[2 2],[1 2])  

if nargin==3
  offset=0;
end

if nargin<=4
  lims=[Sc(1)*Sd(1)*Sn(1) Sc(2)*Sd(2)*Sn(2)];
end

if lims(1)>Sc(1)*Sd(1)*Sn(1)
 error('lims(1) must be <= Sc(1)*Sd(1)*Sn(1)')
end

drawgrid(Sc,Sd,Sn)

ncnt=-1;
if offset==0
  str='';
else
  str=sprintf('%i',0);
  for i=1:offset-1
    str=[str ',' sprintf('%i',i)];
  end
end

for jn=0:Sn(2)-1
 for in=0:Sn(1)-1
  for jd=0:Sd(2)-1
   for id=0:Sd(1)-1
    for jc=0:Sc(2)-1
     if jc+Sc(2)*(jd+Sd(2)*jn)<lims(2)
      for ic=0:Sc(1)-1
       if ic+Sc(1)*(id+Sd(1)*in)<lims(1)
        ncnt=ncnt+1;
        % Physical core location
        i=ijk2i(Sc(1),Sd(1),ic,id,in);
        j=ijk2i(Sc(2),Sd(2),jc,jd,jn);
        % Geographic rank assuming no gaps f(i,j)
      % gpe=i+Sc(1)*Sd(1)*Sn(1)*j;
        gpe=i+lims(1)*j;
        % Geographic (i,j) based on used core count
      % [ii,jj]=i2ij(Sc(1)*Sd(1)*Sn(1),ncnt);
      % disp([i ii j jj ii-i jj-j])
      % disp([jn in jd id jc ic ncnt i j gpe])
        text(i,j,sprintf('%i',gpe+offset),'HorizontalAlignment','Center')
        text(i+.25,j-.25,sprintf('(%i)',ncnt+offset),'HorizontalAlignment','Center','FontSize',6)
        if length(str)==0
          str=sprintf('%i',gpe+offset);
        else
          str=[str ',' sprintf('%i',gpe+offset)];
        end
       end % if lims(1)
      end % ic
     end % if lims(2)
    end % jc
   end % id
  end % jd
 end % in
end % jn
tag=sprintf('%ix%i,%ix%i,%ix%i=%ix%i=%i',[Sc Sd Sn],Sc.*Sd.*Sn,prod(Sc.*Sd.*Sn));
file=['MPICH_RANK_ORDER.' tag];
disp(tag)
disp(str)
if strcmp(input('Write MPICH_RANK_ORDER (y/n)?\n','s'),'y')
  disp( sprintf('Writing %s ...',file) )
  fid=fopen(file,'wt');
  fprintf(fid,'%s\n',str)
  fclose(fid)
  disp('...done!')
end

% =============================================================================

function [i] = ijk2i(C,D,ic,id,in)
i=ic+C*(id+D*in);

% =============================================================================

function [i,j] = i2ij(C,n)
i=mod(n,C);
j=round((n-i)/C); % round should be unnecessary

% =============================================================================

function [] = drawgrid(Sc,Sd,Sn)
clf
nx=Sc(1)*Sd(1)*Sn(1);
ny=Sc(2)*Sd(2)*Sn(2);

for j=0:ny
 h=line([0 nx]-.5,[1 1]*j-.5);
 set(h,'Color',[0 0 0])
end
for i=0:nx
 h=line([1 1]*i-.5,[0 ny]-.5);
 set(h,'Color',[0 0 0])
end
axis([0 nx 0 ny]-.5)
set(gca,'XTick',0:nx,'YTick',0:ny)

for j=0:Sc(2):ny
 h=line([0 nx]-.5,[1 1]*j-.5);
 set(h,'Color',[1 0 0],'LineWidth',2)
end
for i=0:Sc(1):nx
 h=line([1 1]*i-.5,[0 ny]-.5);
 set(h,'Color',[1 0 0],'LineWidth',2)
end

for j=0:Sc(2)*Sd(2):ny
 h=line([0 nx]-.5,[1 1]*j-.5);
 set(h,'Color',[0 0 1],'LineWidth',2)
end
for i=0:Sc(1)*Sd(1):nx
 h=line([1 1]*i-.5,[0 ny]-.5);
 set(h,'Color',[0 0 1],'LineWidth',2)
end
