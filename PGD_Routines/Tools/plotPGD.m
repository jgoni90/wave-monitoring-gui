function plotPGD(pgd,option,colorspec,step,maxTerms)

if ~iscell(pgd.errors.residual)
    aux = pgd.errors.residual;
    pgd.errors.residual = cell(1,1);
    pgd.errors.residual{1} = aux;
end

if nargin == 3
    step = 1;
    maxTerms = size(pgd.RB{1},2);
elseif nargin == 4
    maxTerms = size(pgd.RB{1},2);
end

nOfProjections = length(pgd.errors.residual);

fin = 1;
ivar0 = [];
for i = 1:nOfProjections
    iterms = length(pgd.errors.residual{i});
    if strcmp(option,'BYPROJECTION') && i > 1
        iprojectedTerms = length(pgd.projection.errors.residual{i-1});
        ivar0 = pgd.projection.errors.residual{i-1}(end);
        ini = iprojectedTerms;
        fin = ini + iterms;
    elseif strcmp(option,'') || i == 1
        ini = fin;
        fin = fin + iterms;
        if i == 1, fin = fin - 1; end
    end
    
    if fin > maxTerms, fin = maxTerms; end
    ivarY = pgd.errors.residual{i}(1:step:end);
    ivarY = [ivar0 ; ivarY];
    ivarX = ini:step:fin;
    ivarY = ivarY(1:length(ivarX));
    posPOS = ivarY > 0;
    ivarYPOS = ivarY(posPOS);
    ivarXPOS = ivarX(posPOS);
    h = semilogy(ivarXPOS,ivarYPOS,colorspec);
    hold on
    
    ivar0 = ivarYPOS(end);
    set(h,'markersize',20,'linewidth',3)
    if fin > maxTerms, break, end
end

set(gca,'fontsize',20)

    
        
      
