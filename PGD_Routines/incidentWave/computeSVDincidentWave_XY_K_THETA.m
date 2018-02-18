function ui = computeSVDincidentWave_XY_K_THETA(PGDmeshes,algorithm,parameters)

%% INCIDENT WAVE MATRICES and MESHES
disp('  Building IW matrices and meshes...')

%Dimensions
dimXY           = findPGDdimension('XY',parameters.PGDdimensions);
dimK            = findPGDdimension('K',parameters.PGDdimensions);
dimTHETA        = findPGDdimension('THETA',parameters.PGDdimensions);
uidimXY         = dimXY;
uidimKTHETA     = dimK;

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

%% COMPUTATION OF SVD

%Parameters
nOfSVDfunctions = 2;
tolSVD          = parameters.iwparam.relerror;
nOfNodesK       = size(PGDmeshes(dimK).X,1);
nOfNodesTHETA   = size(PGDmeshes(dimTHETA).X,1);
nOfNodesXY      = numel(ui.meshes(uidimXY).nodes);
bottom          = parameters.meshes(dimXY).bottom(ui.meshes(uidimXY).nodes(1));
K = zeros(nOfNodesK,1); for i = 1:nOfNodesK, K(i) = computeWaveNumber(PGDmeshes(dimK).X(i),bottom); end

%SVD data
disp('  Computation of matrix data...')

SVDdata = zeros(nOfNodesXY,nOfNodesK,nOfNodesTHETA);
for n = 1:nOfNodesK
    kn = K(n);
    for m = 1:nOfNodesTHETA
        thetam = PGDmeshes(dimTHETA).X(m);
        SVDdata(:,n,m) = exp(1i * kn * (ui.meshes(uidimXY).X(:,1)*cos(thetam) +...
            ui.meshes(uidimXY).X(:,2)*sin(thetam)));
    end
end
SVDdata = reshape(SVDdata,nOfNodesXY,nOfNodesK*nOfNodesTHETA);

%SVD calculation
disp('  Computation of SVD...')

[U,s,V] = svd(SVDdata,0); %economic size SVD
sdiag = diag(s);
nOfTerms = find(sdiag/norm(sdiag) < tolSVD,1);
vterms = 1:nOfTerms;

%Storage (only necessary terms)
disp(['  ' num2str(nOfTerms) ' terms for an alpha tolerance ' num2str(tolSVD)])
ui.nOfTerms = nOfTerms;
ui.alpha = sdiag(vterms).';
ui.RB = cell(nOfSVDfunctions,1);
ui.RB{dimXY} = ui.alpha(ones(nOfNodesXY,1),:) .* U(:,vterms);
ui.RB{dimK} = conj(V(:,vterms));

        
        
        
        
        