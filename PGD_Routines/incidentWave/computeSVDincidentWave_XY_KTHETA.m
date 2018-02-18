function ui = computeSVDincidentWave_XY_KTHETA(PGDmeshes,algorithm,parameters)

%% INCIDENT WAVE MATRICES and MESHES
disp('  Building IW matrices and meshes...')

%Dimensions
dimXY           = findPGDdimension('XY',parameters.PGDdimensions);
dimKTHETA       = findPGDdimension('KTHETA',parameters.PGDdimensions);
uidimXY         = dimXY;
uidimKTHETA     = dimKTHETA;

%Meshes
ui.meshes = struct();
[ui.meshes(uidimXY).X,ui.meshes(uidimXY).T,ui.meshes(uidimXY).nodes] = ...
    getSubmesh(PGDmeshes(dimXY).X,PGDmeshes(dimXY).T.ext);
ui.meshes(uidimXY).referenceElement = PGDmeshes(dimXY).referenceElement;

ui.meshes(uidimKTHETA).X = PGDmeshes(dimKTHETA).X;
ui.meshes(uidimKTHETA).T = PGDmeshes(dimKTHETA).T;
ui.meshes(uidimKTHETA).referenceElement = PGDmeshes(dimKTHETA).referenceElement;

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
nOfNodesKTHETA  = size(PGDmeshes(dimKTHETA).X,1);
nOfNodesXY      = numel(ui.meshes(uidimXY).nodes);
bottom          = parameters.meshes(dimXY).bottom(ui.meshes(uidimXY).nodes(1));

%SVD data
disp('  Computation of matrix data...')

K = zeros(parameters.meshes(dimKTHETA).nK+1,1);
for n = 1:parameters.meshes(dimKTHETA).nK+1
    K(n) = computeWaveNumber(ui.meshes(dimKTHETA).X(n,1),bottom);
end
K = K(:,ones(1,parameters.meshes(dimKTHETA).nTHETA+1));
COEF1 = K(:) .* cos(ui.meshes(dimKTHETA).X(:,2));
COEF2 = K(:) .* sin(ui.meshes(dimKTHETA).X(:,2));
SVDdata = exp(1i * (ui.meshes(uidimXY).X(:,1) * COEF1.' + ui.meshes(uidimXY).X(:,2) * COEF2.'));

%SVD calculation
disp('  Computation of SVD...')

[U,s,V] = svd(SVDdata,0); %economic size SVD
sdiag = diag(s);
nOfTerms = find(sdiag/norm(sdiag) < tolSVD,1);
if isempty(nOfTerms), nOfTerms = numel(sdiag); end
vterms = 1:nOfTerms;

%Storage (only necessary terms)
disp(['  ' num2str(nOfTerms) ' terms for an alpha tolerance ' num2str(tolSVD)])
ui.nOfTerms = nOfTerms;
ui.alpha = sdiag(vterms).';
ui.RB = cell(nOfSVDfunctions,1);
ui.RB{dimXY} = ui.alpha(ones(nOfNodesXY,1),:) .* U(:,vterms);
ui.RB{dimKTHETA} = conj(V(:,vterms));

        
        
        
        
        