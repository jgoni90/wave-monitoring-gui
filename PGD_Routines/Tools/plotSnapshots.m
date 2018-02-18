
%% OPTIONS

interpDim             = {'K' 'THETA'};
paramValues           = {PGDmeshes(2).X(1) PGDmeshes(3).X(1)};
% paramValues          = [0.61,3.3999]; %BCN T1
% paramValues          = [0.42,3.3999]; %BCN T2
% paramValues          = [0.55,3.9]; %BCN T3
% paramValues          = [0.61,3.3929]; %Mataro
runFEM                = 1;
nOfPGDterms           = [];
PGDfile               = '';
PGDdatafile           = '';
GUIfile               = '';

%% DATA

if ~isempty(PGDfile),       load(PGDfile),      end
if ~isempty(PGDdatafile),   load(PGDdatafile),  end
if ~isempty(GUIfile),       load(GUIfile),      end

%% PARAMETER VALUES

dims = findPGDdimension(interpDim,parameters.PGDdimensions);
dimXY = findPGDdimension('XY',parameters.PGDdimensions);
snapvalue = [];

%Freq. and incident angle
if mystrcmp({'K' 'THETA'},interpDim)
    snapvalue(1) = paramValues{1};
    snapvalue(2) = paramValues{2};
    
elseif mystrcmp('K',interpDim)
    snapvalue(1) = paramValues{1};
    snapvalue(2) = parameters.fixedParametricDims(2);
    
elseif mystrcmp('THETA',interpDim)
    snapvalue(1) = parameters.fixedParametricDims(1);
    snapvalue(2) = paramValues{1};
    
elseif mystrcmp('KTHETA',interpDim)
    snapvalue(1) = paramValues{1}(1);
    snapvalue(2) = paramValues{1}(2);
else
    snapvalue = parameters.fixedParametricDims(1:2);
end

%Absorbing coefficients
dimALPHA = findPGDdimension('ALPHA*',parameters.PGDdimensions);
snapALPHA = parameters.fixedParametricDims(3:end);
cont = 1;
for i = intersect(dims,dimALPHA)
    p = find(strcmp(parameters.PGDdimensions{i,1},interpDim));
    snapALPHA(cont) = paramValues{p};
    cont = cont + 1;
end
snapvalue = [snapvalue snapALPHA];

%% COMPUTE AND PLOT SNAPSHOT

%PGD snapshot
if isempty(nOfPGDterms), nOfPGDterms = size(pgd.RB{1},2); end
u = interpolatePGD(pgd,dims,paramValues,PGDmeshes,nOfPGDterms);
    
%FEM snapshot
if runFEM
    ufem = computeFEMfromPGDdata(snapvalue,PGDmeshes,parameters);
    figure
%     plotSolution(PGDmeshes(dimXY).X,PGDmeshes(dimXY).T.int,abs(ufem),PGDmeshes(dimXY).referenceElement);
    plotSolutionInterpElem(data.mesh,abs(ufem),'FEM',20,'intT');
    set(gca,'fontsize',20)
    title('FEM')
    axis off
    acaxis = caxis(gca);
end

%Plot PGD solution
figure
% plotSolution(PGDmeshes(dimXY).X,PGDmeshes(dimXY).T.int,abs(u),PGDmeshes(dimXY).referenceElement);
plotSolutionInterpElem(data.mesh,abs(u),'FEM',20,'intT');
set(gca,'fontsize',20)
title('PGD')
axis off
if runFEM || exist('acaxis','var'), caxis(acaxis), end

