
%% General

clear all
home
parameters.nOfPGDdimensions = 2;
nOfSeriesTerms = 6;

%% 1D dimensions

lim = {[1.5,2],[1.5,2]};
ininOfNodes = 100;
nDeg = 1;

%% Data for PGD projection

parameters.maxterms                 = 500;
algorithm.projection.residualType   = 'ALPHACOEF';
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
pgdtensor.counters.Uterm = nOfSeriesTerms - 1;
[XX,YY] = meshgrid(meshes(1).X,meshes(2).X);
tensor = ones(size(meshes(1).X,1),size(meshes(2).X,1));
for i = 1:pgdtensor.counters.Uterm
    tensor = tensor + (XX - YY).^i;
end
surf(XX,YY,tensor)

%PGD projection
printPGD()
pgd = PGDprojection_NonSeparableTensor(...
    tensor,...
    matrices,...
    meshes,...
    tensor,...
    parameters,...
    algorithm);








