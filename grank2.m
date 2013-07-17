function [list] = grank2(nCores,Grid,nPEs)
% grank2  List possible decompositions for given PE count
%
% grank2(nCores,Grid,nPEs)
%
% Lists the possible decompositions for a given core-count, model
% 2d-grid and PE count
%
%  nCores - (scalar) number of cores per node/blade
%  Grid   - (2-vector) number of points in model grid
%  nPEs   - number of PEs to use
%
% e.g. 24-core, 360x210 grid, 72 PEs
% grank2(24,[360 210],72)

if mod(nPEs,nCores) ~= 0
 error('nPEs must be a multiple of nCores. Do help grank2 for more information')
end

titl=0;
list=[];
coreCombinations=myfactor_products(nCores);
npe=nPEs;
bladeCombinations=myfactor_products(npe/nCores);
for b=1:size(bladeCombinations,1)
  for c=1:size(coreCombinations,1)
    x=Grid(1)/(bladeCombinations(b,1)*coreCombinations(c,1));
    y=Grid(2)/(bladeCombinations(b,2)*coreCombinations(c,2));
    if round(x)==x & round(y)==y
      ncxy=bladeCombinations(b,:).*coreCombinations(c,:);
      list(end+1,:)=[npe bladeCombinations(b,:) coreCombinations(c,:) x y ncxy];
      if nargout==0
        if titl==0
         disp(sprintf('PEs\t(Blades)\t(Cores/blade)\t(Points/core)\tCores'))
         titl=1;
        end
        disp(sprintf('%i\t(%i,%i)\t\t(%i,%i)\t\t(%i,%i)\t\t(%i,%i)', list(end,:) ))
      end
    end
  end % c
end % b
if nargout==0
 clear list
end

% ==============================================================================

function [listOfProducts] = myfactor_products(A)

if isprime(A) || A==1
  listOfProducts=A;
  compliment=1;
else
  listOfProducts=local_factor_products(1,factor(A));
  listOfProducts=unique(listOfProducts);
  compliment=A./listOfProducts;
end

listOfProducts=[listOfProducts' compliment'];

% ------------------------------------------------------------------------------

function [listOfProducts] = local_factor_products(listOfProducts, factors)

newProd=prod(factors);
if newProd>1
  listOfProducts=[listOfProducts prod(factors)];
end
for j=1:length(factors)
 facts=factors; facts(j)=[];
 listOfProducts=local_factor_products(listOfProducts, facts);
end
