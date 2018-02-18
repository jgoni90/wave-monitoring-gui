
%% General

clear all
home
parameters.PGDdimensions    = {'T' ; 'A' ; 'B' ; 'C' ; 'D' ; 'E'};
parameters.nOfPGDdimensions = 6;
computePGDonly              = true;

%% 1D dimensions

lim = [0 2*pi
       -1 1
       -3 3
       0 4
       0 5
       1 12];
ininOfNodes = [100,20,20,20,6,20];
nDeg = [1,1,1,1,1,1];

%% Data for PGD projection

parameters.maxterms                 = 300;
algorithm.projection.residualType   = 'ALPHACOEF';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 5;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-10;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = '';

%% 1D meshes

meshes = struct();
nOfNodes = zeros(size(ininOfNodes));
for i = 1:parameters.nOfPGDdimensions
    [meshes(i).X,meshes(i).T] = createPGDMesh('THETA',...
    struct('THETAini',lim(i,1),'THETAend',lim(i,2),'nDeg',nDeg(i),'ininOfNodes',ininOfNodes(i)));
    nOfNodes(i) = size(meshes(i).X,1);
    meshes(i).referenceElement = createReferenceElement(1,(nDeg(i)+1)*(nDeg(i)+2)/2);
end
i = 5; %Dimension D has to include only entire numbers
meshes(i).X = linspace(lim(i,1),lim(i,2),nOfNodes(i))';

%% Mass matrices and vectors

matrices = struct();
for i = 1:parameters.nOfPGDdimensions
    auxones = {ones(size(meshes(i).X,1),1)};
    [Mw,fw] = PGDmassMatrixAndVector1D(...
        meshes(i).X,...
        meshes(i).T,...
        meshes(i).referenceElement,...
        {auxones},...   
        1,...
        ones(size(meshes(i).T,1),1));
    matrices(i).M = Mw{1}{1};
    matrices(i).f = fw{1}(:,1);
end

%% PGD projections

%Data tensor
T = meshes(1).X;
A = reshape(meshes(2).X,1,nOfNodes(2));
B = reshape(meshes(3).X,1,1,nOfNodes(3));
C = reshape(meshes(4).X,1,1,1,nOfNodes(4));
D = reshape(meshes(5).X,1,1,1,1,nOfNodes(5));
E = reshape(meshes(6).X,1,1,1,1,1,nOfNodes(6));
coef1 = bsxfun(@times,A,exp(cos(T)));
coef2 = bsxfun(@times,B,cos(bsxfun(@times,C,T)));
coef3 = bsxfun(@power,sin(bsxfun(@rdivide,T,E)),D);
tensor = bsxfun(@plus,coef1,bsxfun(@plus,-coef2,coef3));

%PGD projection
fhandle = {str2func('butterfly_T'),...
           str2func('butterfly_A'),...
           str2func('butterfly_B'),...
           str2func('butterfly_C'),...
           str2func('butterfly_D'),...
           str2func('butterfly_E')};
printPGD()
pgd = PGDprojection_NonSeparableVector(...
    fhandle,...
    matrices,...
    meshes,...
    tensor,...
    parameters,...
    algorithm);

if computePGDonly, return, end

%% Execute for interpolation and plot results

%Interpolation nodes for parameters
pos = [7 4 1 3 19]; 
%pos = [20 15 20 6 20];
%pos = [10 2 10 4 12]; 

%PGD solution r(THETA)
vterms = 1:274;
fi = pgd.RB{2}(pos(1),vterms).*pgd.RB{3}(pos(2),vterms).*pgd.RB{4}(pos(3),vterms).*...
     pgd.RB{5}(pos(4),vterms).*pgd.RB{6}(pos(5),vterms);
r = pgd.RB{1}(:,vterms)*fi.';

%Exact solution r_(THETA)
A_ = A(pos(1)); B_ = B(pos(2)); 
C_ = C(pos(3)); D_ = D(pos(4)); 
E_ = E(pos(5));
r_ = A_.*exp(cos(T)) - B_.*cos(C_.*T) + (sin(T./E_)).^D_;

%Plot
polar(T,r,'b-')
hold on
polar(T,r_,'ro')

%% PGD error in the tensor

epgd = zeros(pgd.counters.Uterm,1);
pgdtensor = 0;
for i = 1:pgd.counters.Uterm
    Tp = pgd.RB{1}(:,i);
    Ap = reshape(pgd.RB{2}(:,i),1,nOfNodes(2));
    Bp = reshape(pgd.RB{3}(:,i),1,1,nOfNodes(3));
    Cp = reshape(pgd.RB{4}(:,i),1,1,1,nOfNodes(4));
    Dp = reshape(pgd.RB{5}(:,i),1,1,1,1,nOfNodes(5));
    Ep = reshape(pgd.RB{6}(:,i),1,1,1,1,1,nOfNodes(6));
    pgdtensor = pgdtensor + ...
        bsxfun(@times,bsxfun(@times,bsxfun(@times,bsxfun(@times,bsxfun(@times,Ep,Dp),Cp),Bp),Ap),Tp);
    epgd(i) = norm(pgdtensor(:) - tensor(:)) / norm(tensor(:));
    disp(epgd(i))
end

%% HOSVD model

%Number of terms of the tucker decomposition for each case and dimension
hosvdterms = [
              1 1 1 1 1 1
              2 2 2 2 2 2
              3 3 3 3 3 3
              4 4 4 4 4 4
              5 5 5 5 5 5
              6 6 6 6 6 6
              10 10 10 10 6 10
              15 15 15 15 6 15
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









