function ui = computePGDincidentWave_XY_K_THETA(PGDmeshes,algorithm,parameters)

%% INCIDENT WAVE MATRICES and MESHES
disp('  Building IW matrices and meshes...')

%Dimensions
dimXY           = findPGDdimension('XY',parameters.PGDdimensions);
dimK            = findPGDdimension('K',parameters.PGDdimensions);
dimTHETA        = findPGDdimension('THETA',parameters.PGDdimensions);
uidimXY         = dimXY;
uidimKTHETA     = dimK;
nOfuidimensions = 2;

%Meshes
ui.meshes = struct();
[ui.meshes(uidimXY).X,ui.meshes(uidimXY).T,ui.meshes(uidimXY).nodes] = ...
    getSubmesh(PGDmeshes(dimXY).X,PGDmeshes(dimXY).T.ext);
ui.meshes(uidimXY).referenceElement = PGDmeshes(dimXY).referenceElement;
[ui.meshes(uidimKTHETA).X,ui.meshes(uidimKTHETA).T] = createPGDMesh('DoQUAD',...
        struct('x',PGDmeshes(dimK).X,'y',PGDmeshes(dimTHETA).X));
ui.meshes(uidimKTHETA).referenceElement = createReferenceElement(0,4,[]);

%Mass matrices
ui.matrices = struct();
for i = [uidimXY,uidimKTHETA]
    auxones = ones(size(ui.meshes(i).X,1),1);
    [~,Mint] = PGDberkhoffVolumeMatrices(...
        ui.meshes(i).X,...
        ui.meshes(i).T,...
        ui.meshes(i).referenceElement,...
        auxones,...
        auxones,...
        auxones,...
        auxones,...
        []); %no PML elements
    ui.matrices(i).M = Mint{1};
end

%% DATA
disp('  Dimensionless IW in X and Y directions...')

%Load dimensionless 1D pgd for X and Y directions
pgdX = load(parameters.iwparam.loadfileX); 
pgdY = load(parameters.iwparam.loadfileY);

%Some variables
algorithmcopy   = algorithm;
parameterscopy  = parameters;
xnodes          = PGDmeshes(dimXY).X(ui.meshes(uidimXY).nodes,1);
xmax            = max(xnodes);
xmin            = min(xnodes);
Dx              = xmax - xmin;
ynodes          = PGDmeshes(dimXY).X(ui.meshes(uidimXY).nodes,2);
ymax            = max(ynodes);
ymin            = min(ynodes);
Dy              = ymax - ymin;
bottom          = parameters.meshes(dimXY).bottom(ui.meshes(uidimXY).nodes(1));

%Function handles
oneDimMap               = @(d1,d2,x)(1/(d1(1)-d1(2)))*((d2(1)-d2(2))*x ...
                                    + d2(2)*d1(1) - d2(1)*d1(2));
omegaFromKdimensionless = @(k)sqrt(k.*9.81.*tanh(k));

%Dimensionless wavenumber and frequency
nOfNodesK = size(PGDmeshes(dimK).X,1);
KDx       = zeros(nOfNodesK,1);
KDy       = zeros(nOfNodesK,1);
for i = 1:nOfNodesK
    k = computeWaveNumber(PGDmeshes(dimK).X(i),bottom);
    KDx(i) = k * Dx;
    KDy(i) = k * Dy; 
end
freqx = omegaFromKdimensionless(KDx);
freqy = omegaFromKdimensionless(KDy);

%% INTERPOLATION on the PML area
disp('  Interpolation on PGD meshes...')

%Interpolate PGD of the incident wave in X x (OMEGA x THETA)
uiX.RB = cell(pgdX.parameters.nOfPGDdimensions,1);
uiX.counters = pgdX.pgdproj.counters;

uiX.RB{uidimXY} = interpolateOn1Dmesh(pgdX.pgdproj.RB{uidimXY},xnodes,...
    oneDimMap([0,1],[xmin,xmax],pgdX.meshes(uidimXY).X),...
    pgdX.meshes(uidimXY).T,pgdX.meshes(uidimXY).referenceElement);

uiX.RB{uidimKTHETA} = projectQUADonQUADmesh(pgdX.pgdproj.RB{uidimKTHETA},...
    {freqx,PGDmeshes(dimTHETA).X},pgdX.meshes(uidimKTHETA).X,...
    pgdX.meshes(uidimKTHETA).T,pgdX.meshes(uidimKTHETA).referenceElement);

clear pgdX

%Interpolate PGD of the incident wave in Y x (OMEGA x THETA)
uiY.RB = cell(pgdY.parameters.nOfPGDdimensions,1);
uiY.counters = pgdY.pgdproj.counters;

uiY.RB{uidimXY} = interpolateOn1Dmesh(pgdY.pgdproj.RB{uidimXY},ynodes,...
    oneDimMap([0,1],[ymin,ymax],pgdY.meshes(uidimXY).X),...
    pgdY.meshes(uidimXY).T,pgdY.meshes(uidimXY).referenceElement);

uiY.RB{uidimKTHETA} = projectQUADonQUADmesh(pgdY.pgdproj.RB{uidimKTHETA},...
    {freqy,PGDmeshes(dimTHETA).X},pgdY.meshes(uidimKTHETA).X,...
    pgdY.meshes(uidimKTHETA).T,pgdY.meshes(uidimKTHETA).referenceElement.degree);

clear pdgY

%% PGD PROJECTION
disp('  PGD projection in XY x (OMEGA x THETA)...')

%Projection parameters
parameterscopy.nOfPGDdimensions         = nOfuidimensions;
parameterscopy.maxterms                 = parameters.iwparam.maxterms;
algorithmcopy.projection.residualType   = parameters.iwparam.resisualType;
parameterscopy.residualEachTerm         = 1;
algorithmcopy.projection.maxiter        = parameters.iwparam.maxiter;
parameterscopy.toliter                  = 1e-3;
algorithmcopy.projection.relerror       = parameters.iwparam.relerror;
algorithmcopy.projection.projCoeff      = true;

%Do projection in XY x (OMEGA x THETA)
printPGD()
pgdproj = PGDprojection_separableFIRSTsingleProduct...
            ({uiY,uiX},...
            ui.matrices,...
            ui.meshes,...
            [],...
            parameterscopy,...
            algorithmcopy,...
            parameters.iwparam.maxGBmemory);
clear uiY uiX
ui.RB = pgdproj.RB;

        
        
        
        
        