home
clear all
close all

%% SET DATA
pgdfile = 'RB_mataro4D_i3_pg.mat';
datafile = 'PGDdata_mataro4D.mat';

interpPlot = true;
interpDeg = 20;
dataGUI = 'data_P08.mat';

refname = 'RB_mataro4D_i3_pg_VIDEOIMAG';

fixedPointArrow = [-750,-250];
lengthArrow = 200;
rotPoint = [-300,-270];

duration = 7;
dt = 0.1; %minimum in iMovie

%% PARAMETRIC DATA
load(pgdfile)
load(datafile)
if interpPlot, load(dataGUI), end

nOfTerms = size(pgd.RB{1},2);
nOfImages = round(duration/dt);

paramDims = findPGDdimension('PARAMETRICDIM',parameters.PGDdimensions);
spaceDim = setdiff(1:parameters.nOfPGDdimensions,paramDims);
nOfParamDims = numel(paramDims);
interpX = zeros(nOfParamDims,nOfImages);
indexParam = zeros(1,parameters.nOfPGDdimensions);
for i = 1:nOfParamDims
    iparam = paramDims(i);
    ini = PGDmeshes(iparam).X(1);
    fin = PGDmeshes(iparam).X(end);
    interpX(i,:) = linspace(ini,fin,nOfImages);
    indexParam(iparam) = i;
end

%% INTERPOLATED SPACIAL SOLUTION
u = zeros(size(PGDmeshes(spaceDim).X,1),nOfImages);
for frame = 1:nOfImages
    if nOfParamDims == 2
        u(:,frame) = interpolatePGD(2+nOfParamDims,...
                interpX(1,frame),interpX(2,frame),...
                pgd.RB{paramDims(1)},pgd.RB{paramDims(2)},pgd.RB{spaceDim},...
                PGDmeshes(paramDims(1)).X,PGDmeshes(paramDims(1)).T,parameters.meshes(paramDims(1)).nDeg,...
                PGDmeshes(paramDims(2)).X,PGDmeshes(paramDims(2)).T,parameters.meshes(paramDims(2)).nDeg,...
                ones(1,nOfTerms));
    elseif nOfParamDims == 1
        u(:,frame) = interpolatePGD(2+nOfParamDims,...
                interpX(1,frame),...
                pgd.RB{paramDims(1)},pgd.RB{spaceDim},...
                PGDmeshes(paramDims(1)).X,PGDmeshes(paramDims(1)).T,parameters.meshes(paramDims(1)).nDeg,...
                ones(1,nOfTerms));
    end
    %%%%%%% FOR PLOTTING THE WAVE AMPLIFICATION FACTOR (comment/uncomment)
    kvector = interpX(1,frame)*[cos(interpX(2,frame)) sin(interpX(2,frame))];
    u0 = exp(sqrt(-1)*(kvector(1)*data.mesh.X(:,1) + kvector(2)*data.mesh.X(:,2)));
    u(:,frame) = abs(u(:,frame) + u0);
    %%%%%%%
end
maxval = max(max(real(u)));
minval = min(min(real(u)));
axval = [minval maxval];

%% CREATE AND SAVE IMAGES
thetaDim = findPGDdimension('THETA',parameters.PGDdimensions);
for frame = 1:nOfImages
    
    framename = [refname num2str(frame)];
    
    %Plot
    if interpPlot
        H = plotSolutionInterpElem(data.mesh,u(:,frame),'FEM',interpDeg);
    else
        H = plotSolution(PGDmeshes(spaceDim).X,PGDmeshes(spaceDim).T,u(:,frame),PGDmeshes(spaceDim).referenceElement);
    end
    axis off
    caxis(axval)
    delete(H(end))
    
    %Arrow
    hold on
    Ha = arrow('start',[fixedPointArrow(1)-lengthArrow fixedPointArrow(2)],'stop',fixedPointArrow,'baseangle',90);
    if thetaDim
        angle = interpX(indexParam(thetaDim),frame);
    else
        angle = parameters.fixedParametricDims(2);
    end
    rotate(Ha,[0 0 1],180*angle/pi,[rotPoint 0]);

    %Image
    print([framename '.jpg'],'-djpeg')
    
    delete([H(1:end-1) Ha])
end
            
    