
%% General

clear all
home
parameters.nOfPGDdimensions = 2;
nOfSeriesTerms = 10;

%% 1D dimensions

lim = {[0,1],[0,1]};
x0 = 0;
y0 = 0;
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
[XX,YY] = meshgrid(meshes(1).X,meshes(2).X);
pgdtensor.counters.Uterm = nOfSeriesTerms^2;
tensor = 0;
iterm = 1;
polx = 0;
for i = 0:nOfSeriesTerms-1
    coef = (meshes(1).X - x0).^i;
    polx = polx + coef; %Polynomial vector Px, assuming Px = Py
    for j = 0:nOfSeriesTerms-1
        pgdtensor.RB{1}(:,iterm) = coef;
        pgdtensor.RB{2}(:,iterm) = (meshes(2).X - y0).^j;
        tensor = tensor + pgdtensor.RB{1}(:,iterm) * (pgdtensor.RB{2}(:,iterm)).';
        iterm = iterm + 1;
    end
end

%Plot tensor
h = surfc(XX,YY,tensor,'edgealpha',0,'facealpha',0.5);

%PGD projection
printPGD()
pgd = PGDprojection(pgdtensor,...
    matrices,...
    meshes,...
    tensor,...
    parameters,...
    algorithm);








