
%% General

clear all
home
parameters.PGDdimensions    = {'X' ; 'K'};
parameters.nOfPGDdimensions = 2;

%% Space / wave number / angle data

direction       = 1; %1 --> exp(i k x cos(theta)) and 2 --> exp(i k y sin(theta))

aX              = 0;
bX              = 1;
aK              = 1;
bK              = 50;
aTH             = 0;
bTH             = pi;
aKp             = -50;
bKp             = 50;

ininOfNodes_x   = 300;
nDeg_x          = 1;
ininOfNodes_w   = 1000;
nDeg_w          = 1;
ininOfNodes_TH  = 1000;
nDeg_TH         = 1;
ininOfNodes_Kp  = 300;
nDeg_Kp         = 1;

%% Data for PGD projection

parameters.maxterms                 = 350;
algorithm.projection.residualType   = 'TENSOR_APP';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 20;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-8;
algorithm.projection.projCoeff      = false;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = 'PGD_UI_OUTPUT.mat';

%% Meshes dimension SPACE, WAVE NUMBER, THETA and KP

[meshes(1).X,meshes(1).T] = createPGDMesh('THETA',...
struct('THETAini',aX,'THETAend',bX,'nDeg',nDeg_x,'ininOfNodes',ininOfNodes_x));
nOfNodes_x = size(meshes(1).X,1);
meshes(1).referenceElement = createReferenceElement(1,(nDeg_x+1)*(nDeg_x+2)/2);

[meshes(2).X,meshes(2).T] = createPGDMesh('THETA',...
    struct('THETAini',aK,'THETAend',bK,'nDeg',nDeg_w,'ininOfNodes',ininOfNodes_w));
nOfNodes_w = size(meshes(2).X,1);
meshes(2).referenceElement = createReferenceElement(1,(nDeg_w+1)*(nDeg_w+2)/2);

[meshes(3).X,meshes(3).T] = createPGDMesh('THETA',...
    struct('THETAini',aTH,'THETAend',bTH,'nDeg',nDeg_TH,'ininOfNodes',ininOfNodes_TH));
nOfNodes_TH = size(meshes(3).X,1);
meshes(3).referenceElement = createReferenceElement(1,(nDeg_TH+1)*(nDeg_TH+2)/2);

[meshes(4).X,meshes(4).T] = createPGDMesh('THETA',...
    struct('THETAini',aKp,'THETAend',bKp,'nDeg',nDeg_Kp,'ininOfNodes',ininOfNodes_Kp));
nOfNodes_Kp = size(meshes(4).X,1);
meshes(4).referenceElement = createReferenceElement(1,(nDeg_Kp+1)*(nDeg_Kp+2)/2);

%% Mass matrices dimension SPACE, WAVE NUMBER, THETA and Kp

for i = 1:4 
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

%% PGD projection to decompose in X and Kp

%Data tensor
t = zeros(nOfNodes_x,nOfNodes_Kp);
for j = 1:nOfNodes_Kp
    kp = meshes(4).X(j);
    t(:,j) = exp(1i*kp*meshes(1).X);
end

%Shape functions tensor for PGD
A = tensor_tijNiNj(meshes(1).X,meshes(1).T,meshes(4).X,meshes(4).T,...
    meshes(1).referenceElement,meshes(4).referenceElement,t);

%PGD projection
printPGD()
pgdKp = PGDprojection_NonSeparableTensor(...
            A,...
            matrices([1,4]),...
            meshes([1,4]),...
            t,...
            parameters,...
            algorithm);
clear t A
        
%% Decompose each function of Kp in K x THETA

THETAfun = {@(x)cos(x), @(y)sin(y)};
j = 120;

%Interpolate data tensor
kptensor = meshes(2).X * THETAfun{direction}(meshes(3).X)';
tG = interpolateOn1Dmesh(pgdKp.RB{2}(:,j),kptensor(:),...
    meshes(4).X,meshes(4).T,meshes(4).referenceElement);
tG = reshape(tG,nOfNodes_w,nOfNodes_TH);

%Shape functions tensor for PGD
AG = tensor_tijNiNj(meshes(2).X,meshes(2).T,meshes(3).X,meshes(3).T,...
    meshes(2).referenceElement,meshes(3).referenceElement,tG);

%PGD projection
printPGD()
pgdG = PGDprojection_NonSeparableTensor(...
            AG,...
            matrices([2,3]),...
            meshes([2,3]),...
            tG,...
            parameters,...
            algorithm);


        












