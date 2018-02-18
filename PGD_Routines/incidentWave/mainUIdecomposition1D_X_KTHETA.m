
%% General

clear all
home
parameters.PGDdimensions    = {'X' ; 'KTHETA'};
parameters.nOfPGDdimensions = 2;

%% Space / wave number / angle data

direction       = 2; %1 --> exp(i k x cos(theta)) and 2 --> exp(i k y sin(theta))

aX              = 0;
bX              = 1;
aK              = 1;
bK              = 600;
aTH             = pi;
bTH             = 2*pi;

ininOfNodes_x   = 500;
nDeg_x          = 1;
ininOfNodes_w   = 100;
ininOfNodes_TH  = 100;
nDeg_w_TH       = 1;

%% Data for PGD projection

parameters.maxterms                 = 350;
algorithm.projection.residualType   = 'ALPHACOEF';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 5;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-8;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = 'PGD_UIX_X_KTHETA_1000_500x500[0to1][1to600][PIto2PI].mat';

%% Meshes dimension SPACE, OMEGA x THETA

[meshes(1).X,meshes(1).T] = createPGDMesh('THETA',...
struct('THETAini',aX,'THETAend',bX,'nDeg',nDeg_x,'ininOfNodes',ininOfNodes_x));
nOfNodes_x = size(meshes(1).X,1);
meshes(1).referenceElement = createReferenceElement(1,(nDeg_x+1)*(nDeg_x+2)/2);

omegaFromKdimensionless = @(k)sqrt(k.*9.81.*tanh(k));
aK = omegaFromKdimensionless(aK);
bK = omegaFromKdimensionless(bK);
[meshes(2).X,meshes(2).T] = createPGDMesh('KTHETA',...
    struct('THETAini',aTH,'THETAend',bTH,'nDeg',nDeg_w_TH,...
    'Kini',aK,'Kend',bK,'nK',ininOfNodes_w-1,'nTHETA',ininOfNodes_TH-1,...
    'equalSpaced',true));
nOfNodes_w_TH = size(meshes(2).X,1);
meshes(2).referenceElement = createReferenceElement(0,(nDeg_w_TH+1)^2);

%% Mass matrices dimension SPACE, WAVE NUMBER x THETA

%Dimension SPACE
auxones = {ones(size(meshes(1).X,1),1)};
    Mw = PGDmassMatrix1D(...
        meshes(1).X,...
        meshes(1).T,...
        meshes(1).referenceElement,...
        {auxones},...
        1,...
        ones(size(meshes(1).T,1),1));
matrices(1).M = Mw{1}{1};

%Dimension OMEGA x THETA
auxones = ones(size(meshes(2).X,1),1);
    [~,M] = PGDberkhoffVolumeMatrices(...
        meshes(2).X,...
        meshes(2).T,...
        meshes(2).referenceElement,...
        auxones,...
        auxones,...
        auxones,...
        auxones,...
        []); %no PML elements!
matrices(2).M = M{1};

%% PGD projection to decompose in X and OMEGA x THETA

%Vector of wave numbers
Knodes = zeros(ininOfNodes_w,ininOfNodes_TH);
for i = 1:ininOfNodes_w, Knodes(i,1) = computeWaveNumber(meshes(2).X(i,1),1); end
Knodes = Knodes(:,ones(1,ininOfNodes_TH));

%Data tensor
THETAfun = {@(x)cos(x), @(y)sin(y)};
kp = Knodes(:) .* THETAfun{direction}(meshes(2).X(:,2));
t = exp(1i * meshes(1).X * kp.');

%Shape functions tensor for PGD
A = tensor_tijkNiNjk(meshes(1).X,meshes(1).T,meshes(2).X,meshes(2).T,...
    meshes(1).referenceElement,meshes(2).referenceElement,t);

%PGD projection
printPGD()
pgdproj = PGDprojection_NonSeparableTensor(...
            A,...
            matrices,...
            meshes,...
            t,...
            parameters,...
            algorithm);
        
%% Save

gb = MSL_getGb(pgdproj);
if gb > 1.95
    save(parameters.outputFileName,'pgdproj','parameters','meshes','matrices','-v7.3')
else
    save(parameters.outputFileName,'pgdproj','parameters','meshes','matrices')
end







