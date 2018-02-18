function ui = computePGDincidentWave_XY_K_THETA2(PGDmeshes,algorithm,parameters)

%Load
pgdX = load('PGD_UIX_100x50x50'); 
pgdY = load('PGD_UIY_100x50x50');
nOfterms = 60;
for i = 1:pgdX.parameters.nOfPGDdimensions
    pgdX.pgdproj.RB{i} = pgdX.pgdproj.RB{i}(:,1:nOfterms);
    pgdY.pgdproj.RB{i} = pgdY.pgdproj.RB{i}(:,1:nOfterms);
end

%Dimensions
dimXY = findPGDdimension('XY',parameters.PGDdimensions);
nOfNodesXY = size(PGDmeshes(dimXY).X,1);
dimK = findPGDdimension('K',parameters.PGDdimensions);
dimTHETA = findPGDdimension('THETA',parameters.PGDdimensions);

%Some variables
algorithmcopy = algorithm;
parameterscopy = parameters;
nodesPML = unique(PGDmeshes(dimXY).T.ext);
xnodes = PGDmeshes(dimXY).X(nodesPML,1);
ynodes = PGDmeshes(dimXY).X(nodesPML,2);

%Function handles
oneDimMap = @(d1,d2,x)(1/(d1(1)-d1(2)))*((d2(1)-d2(2))*x + d2(2)*d1(1) - d2(1)*d1(2));

%PGD of the incident wave in X x OMEGA x K
uiX.RB = cell(parameters.nOfPGDdimensions,1);

uiX.RB{dimXY} = interpolateOn1Dmesh(pgdX.pgdproj.RB{1},xnodes,...
    oneDimMap([0,1],[min(xnodes),max(xnodes)],pgdX.meshes(1).X),...
    pgdX.meshes(1).T,pgdX.meshes(1).referenceElement);

uiX.RB{dimK} = interpolateOn1Dmesh(pgdX.pgdproj.RB{2},PGDmeshes(dimK).X,...
    pgdX.meshes(2).X,pgdX.meshes(2).T,pgdX.meshes(2).referenceElement);

uiX.RB{dimTHETA} = interpolateOn1Dmesh(pgdX.pgdproj.RB{3},PGDmeshes(dimTHETA).X,...
    pgdX.meshes(3).X,pgdX.meshes(3).T,pgdX.meshes(3).referenceElement);

clear pgdX

%PGD of the incident wave in Y x OMEGA x K
uiY.RB = cell(parameters.nOfPGDdimensions,1);

uiY.RB{dimXY} = interpolateOn1Dmesh(pgdY.pgdproj.RB{1},ynodes,...
    oneDimMap([0,1],[min(ynodes),max(ynodes)],pgdY.meshes(1).X),...
    pgdY.meshes(1).T,pgdY.meshes(1).referenceElement);

uiY.RB{dimK} = interpolateOn1Dmesh(pgdY.pgdproj.RB{2},PGDmeshes(dimK).X,...
    pgdY.meshes(2).X,pgdY.meshes(2).T,pgdY.meshes(2).referenceElement);

uiY.RB{dimTHETA} = interpolateOn1Dmesh(pgdY.pgdproj.RB{3},PGDmeshes(dimTHETA).X,...
    pgdY.meshes(3).X,pgdY.meshes(3).T,pgdY.meshes(3).referenceElement);

clear pdgY

%PGD of the incident wave in XY x OMEGA x THETA
f1 = [1,0,0];
f2{1}{1} = nodesPML;
f2{1}{2} = nOfNodesXY;
uiprod = arrangePGD(uiX,uiY,f1,f2,'prod');
uiprod.counters.Uterm = size(uiprod.RB{1},2);
clear uiX uiY

%Mass matrices
auxones = ones(size(PGDmeshes(dimXY).X,1),1);
[~,Mint,~,~,Mpml] = PGDberkhoffVolumeMatrices(...
    PGDmeshes(dimXY).X,...
    PGDmeshes(dimXY).T.all,...
    PGDmeshes(dimXY).referenceElement,...
    auxones,...
    auxones,...
    auxones,...
    auxones,...
    parameters.meshes(dimXY).PML.elements);
matrices(dimXY).M = Mint{1} + Mpml{1};

for i = [dimK,dimTHETA]
    auxones = {ones(size(PGDmeshes(i).X,1),1)};
    Mw = PGDmassMatrix1D(...
        PGDmeshes(i).X,...
        PGDmeshes(i).T,...
        PGDmeshes(i).referenceElement,...
        {auxones},...   
        1,...
        ones(size(PGDmeshes(i).T,1),1));
    matrices(i).M = Mw{1}{1};
end

%PGD projection for reducing the number of terms
parameterscopy.maxterms                 = 200;
algorithmcopy.projection.residualType   = 'ALPHACOEF';
parameterscopy.residualEachTerm         = 1;
algorithmcopy.projection.maxiter        = 20;
parameterscopy.toliter                  = 1e-3;
algorithmcopy.projection.relerror       = 1e-3;
algorithmcopy.projection.projCoeff      = true;

printPGD()
ui = PGDprojection(uiprod,...
            matrices,...
            PGDmeshes,...
            [],...
            parameterscopy,...
            algorithmcopy);
        
        
        
        
        