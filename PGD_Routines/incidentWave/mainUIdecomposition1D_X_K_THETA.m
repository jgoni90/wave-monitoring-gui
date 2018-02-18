
%% General

clear all
home
parameters.PGDdimensions    = {'X' ; 'K' ; 'THETA'};
parameters.nOfPGDdimensions = 3;

%% Space / wave number / angle data

UIstr           = 'Y'; %X --> exp(i k x cos(theta)) and Y --> exp(i k y sin(theta))
UIstrp          = 'ReIm'; %Re --> real(ui)  Im --> imag(ui)  ReIm --> both

aX              = 0;
bX              = 1;
aK              = computeWaveNumber(2*pi/1.01,1);
bK              = computeWaveNumber(2*pi/0.59,1);
aTH             = pi;
bTH             = 3*pi/2;

ininOfNodes_x   = 100;
nDeg_x          = 1;
ininOfNodes_w   = 50;
nDeg_w          = 1;
ininOfNodes_TH  = 50;
nDeg_TH         = 1;

%% Data for PGD projection

parameters.maxterms                 = 200;
algorithm.projection.residualType   = 'TENSOR_APP';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 20;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-3;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = 'PGD_UIY_100x50x50.mat';

%% Meshes dimension SPACE, WAVE NUMBER and THETA

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

%% Mass matrices dimension SPACE, WAVE NUMBER and THETA

for i = 1:3 
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
if strcmpi(algorithm.projection.residualType,'TENSOR_APP')
    if strcmpi(UIstr,'X'), THETAfun = @(x)cos(x); elseif strcmpi(UIstr,'Y'), THETAfun = @(y)sin(y); end 
    t = zeros(nOfNodes_x,nOfNodes_w,nOfNodes_TH);
    for n = 1:nOfNodes_w
        k = meshes(2).X(n);
        for m = 1:nOfNodes_TH
            theta = meshes(3).X(m);
            t(:,n,m) = exp(1i*k*meshes(1).X*THETAfun(theta));
        end
    end
else
    t = [];
end

%PGD projection for real part
if any(strcmpi(UIstrp,{'Re' 'ReIm'}))
    fhandle = {str2func(['projectionReUIVector' UIstr '_NonSeparable_X']),...
               str2func(['projectionReUIVector' UIstr '_NonSeparable_K']),...
               str2func(['projectionReUIVector' UIstr '_NonSeparable_THETA'])};
    printPGD()
    pgdRe = PGDprojection_NonSeparableVector(...
        fhandle,...
        matrices,...
        meshes,...
        real(t),...
        parameters,...
        algorithm);
end

%PGD projection for imag part
if any(strcmpi(UIstrp,{'Im' 'ReIm'}))
    fhandle = {str2func(['projectionImUIVector' UIstr '_NonSeparable_X']),...
               str2func(['projectionImUIVector' UIstr '_NonSeparable_K']),...
               str2func(['projectionImUIVector' UIstr '_NonSeparable_THETA'])};
    printPGD()
    pgdIm = PGDprojection_NonSeparableVector(...
        fhandle,...
        matrices,...
        meshes,...
        imag(t),...
        parameters,...
        algorithm);
end

%PGD structure for complex space
if strcmpi(UIstrp,'ReIm')
    pgd = arrangePGD(pgdRe,pgdIm,[1,1,1],[1i,1,1],'add');

    %PGD projection for complex space using the current separable function pgd
    algorithm.projection.residualType = 'TENSOR_APP';
    algorithm.projection.projCoeff = true;
    algorithm.projection.relerror = 1e-15;
    pgdproj = PGDprojection(pgd,...
                matrices,...
                meshes,...
                t,...
                parameters,...
                algorithm);
end

%% Save

save(parameters.outputFileName,'pgd*','parameters','meshes','matrices')










