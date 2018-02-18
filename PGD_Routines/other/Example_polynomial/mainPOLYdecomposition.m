
%% General

clear all
home
polydeg = 20;
parameters.nOfPGDdimensions = polydeg + 2;

%% 1D dimensions

lim = [0,1];
ininOfNodes = 2;
nDeg = 1;

%% Data for PGD projection

parameters.maxterms                 = 500;
algorithm.projection.residualType   = 'ALPHACOEF';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 15;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-7;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = '';

%% 1D meshes

meshes = struct();
nOfNodes = zeros(1,parameters.nOfPGDdimensions);
for i = 1:parameters.nOfPGDdimensions
    [meshes(i).X,meshes(i).T] = createPGDMesh('THETA',...
    struct('THETAini',lim(1),'THETAend',lim(2),'nDeg',nDeg,'ininOfNodes',ininOfNodes));
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

%Data tensor
data = cell(parameters.nOfPGDdimensions,1);
data{1} = meshes(1).X;
for i = 2:parameters.nOfPGDdimensions
    vpos = [ones(1,i-1),nOfNodes(i)];
    data{i} = reshape(meshes(i).X,vpos);
end
tensor = data{2};
for i = 1:polydeg
    coef1 = bsxfun(@power,data{1},i);
    coef = bsxfun(@times,coef1,data{i+2});
    tensor = bsxfun(@plus,tensor,coef);
end

%PGD tensor
auxones = ones(nOfNodes(1),1);
pgdtensor.counters.Uterm = polydeg + 1;
for i = 1:pgdtensor.counters.Uterm
    pgdtensor.RB{1}(:,i) = meshes(1).X.^(i-1);
end
for j = 2:parameters.nOfPGDdimensions
    pgdtensor.RB{j} = [auxones(:,ones(1,j-2)),meshes(j).X,auxones(:,ones(1,parameters.nOfPGDdimensions-j))];
end

%PGD projection
printPGD()
pgd = PGDprojection(pgdtensor,...
    matrices,...
    meshes,...
    tensor,...
    parameters,...
    algorithm);

return

%% PGD error in the tensor

epgd = zeros(pgd.counters.Uterm,1);
ptensor = 0;
for i = 1:pgd.counters.Uterm
    
    %High-dimensional data of the PGD
    data{1} = pgd.RB{1}(:,i);
    for j = 2:parameters.nOfPGDdimensions
        vpos = [ones(1,j-1),nOfNodes(j)];
        data{j} = reshape(pgd.RB{j}(:,i),vpos);
    end
    coef = data{1};
    for j = 2:parameters.nOfPGDdimensions
        coef = bsxfun(@times,coef,data{j});
    end

    %Error
    ptensor = ptensor + coef;
    epgd(i) = norm(ptensor(:) - tensor(:)) / norm(tensor(:));
    disp(epgd(i))
end

%% HOSVD model

%Number of terms of the tucker decomposition for each case and dimension
hosvdterms = [
              2*ones(1,polydeg+2)
             ];

%Loop in number of test cases
ehosvd = zeros(size(hosvdterms,1),1);
for i = 1:size(hosvdterms,1)

    %Number of factors
    nfac = hosvdterms(i,:);

    %Tucker decomposition
    [factors,core] = tucker(tensor,nfac);

    %Tucker model
    hosvd = nmodel(factors,core);

    %HOSVD error in the tensor
    ehosvd(i) = norm(hosvd(:) - tensor(:)) / norm(tensor(:));
    disp(ehosvd(i))
end









