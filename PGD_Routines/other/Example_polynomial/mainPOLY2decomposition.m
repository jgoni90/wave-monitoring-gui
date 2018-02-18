
%% General

clear all
home
parameters.nOfPGDdimensions = 2;
nOfSeriesTerms = 30;

%% 1D dimensions

lim = {[0,1],[0,1]};
ininOfNodes = 100;
nDeg = 1;

%% Data for PGD projection

parameters.maxterms                 = 500;
algorithm.projection.residualType   = 'TENSOR_APP';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 15;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-15;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = '';

%% 1D meshes

meshes = struct();
nOfNodes = zeros(1,parameters.nOfPGDdimensions);
for i = 1:parameters.nOfPGDdimensions
    [meshes(i).X,meshes(i).T] = createPGDMesh('THETA',...
    struct('THETAini',lim{i}(1),'THETAend',lim{i}(2),'nDeg',nDeg,'ininOfNodes',ininOfNodes));
    nOfNodes(i) = size(meshes(i).X,1);
    meshes(i).referenceElement = createReferenceElement(1,(nDeg+1)*(nDeg+2)/2);
end

%% Mass matrices and vectors

matrices = struct();
for i = 1:parameters.nOfPGDdimensions
    auxones = {ones(size(meshes(i).X,1),1)};
    Mw = PGDmassMatrix1D(...
        meshes(i).X,...
        meshes(i).T,...
        meshes(i).referenceElement,...
        {auxones},...   
        1,...
        ones(size(meshes(i).T,1),1));
    matrices(i).M = Mw{1}{1};
end

%% PGD projections

%PGD tensor
pgdtensor.counters.Uterm = nOfSeriesTerms;
tensor = 0;
for i = 1:pgdtensor.counters.Uterm
    pgdtensor.RB{1}(:,i) = (meshes(1).X).^(i-1);
    pgdtensor.RB{2}(:,i) = (meshes(2).X).^(i-1);
    tensor = tensor + pgdtensor.RB{1}(:,i) * (pgdtensor.RB{2}(:,i)).';
end

%Plot tensor
% [XX,YY] = meshgrid(meshes(1).X,meshes(2).X);
% h = surfc(XX,YY,tensor,'edgealpha',0,'facealpha',0.5);
% set(h,'linewidth',3,'linestyle','-')
% set(gca,'fontsize',30,'LineWidth', 2.5)
% set(gca,'xtick',[0.5,1])
% xlabel('x'), ylabel('y')
% hold on
% [Cplot,hc] = contour3(XX,YY,tensor,3,'k-');
% set(hc,'linewidth',3)

%PGD projection
printPGD()
pgd = PGDprojection(pgdtensor,...
    matrices,...
    meshes,...
    tensor,...
    parameters,...
    algorithm);









