function ui = computeSVDincidentWave_XY_K(PGDmeshes,parameters)

%% INFORMATION OF THE INCIDENT WAVE FUNCTION

dimXY = findPGDdimension('XY',parameters.PGDdimensions);
dimK = findPGDdimension('K',parameters.PGDdimensions);

info.handle         = @uifunction;
info.PGDmeshesField = {'T','ext'};

%% COMPUTATION OF SVD

nOfSVDfunctions = 2;
tolSVD = parameters.iwparam.relerror;
nOfNodes_w = size(PGDmeshes(dimK).X,1);
    
%Spatial nodes
T = myGetField(PGDmeshes(dimXY),info.PGDmeshesField);
xy_nodes = unique(T);
nOfxy_nodes = numel(xy_nodes);
ui.meshes = struct();
ui.meshes(dimXY).nodes = xy_nodes;
ui.meshes(dimK).nodes = [];
bottom = parameters.meshes(dimXY).bottom(xy_nodes(1));

%SVD data
disp('      Computation of matrix data...')
SVDdata = zeros(size(xy_nodes,1),nOfNodes_w);
for j = 1:nOfNodes_w
    SVDdata(:,j) = info.handle(PGDmeshes(dimK).X(j),...
                                   bottom,...
                                   PGDmeshes(dimXY).X(xy_nodes,:),...
                                   parameters.fixedParametricDims(2));
end

%SVD calculation
disp('      Computation of SVD...')
[U,s,V] = svd(SVDdata,0); %economic size SVD
sdiag = diag(s);
nOfTerms = find(sdiag/norm(sdiag) < tolSVD,1);
vterms = 1:nOfTerms;
R = U(:,vterms)*s(vterms,vterms)*V(:,vterms)';
err = norm(SVDdata(:) - R(:)) / norm(R(:));

%Storage (only necessary terms)
disp(['      ' num2str(nOfTerms) ' terms for an alpha tolerance ' num2str(tolSVD) ' with error ' num2str(err)])
ui.nOfTerms = nOfTerms;
ui.alpha = sdiag(vterms).';
ui.RB = cell(nOfSVDfunctions,1);
ui.RB{dimXY} = ui.alpha(ones(nOfxy_nodes,1),:) .* U(:,vterms);
ui.RB{dimK} = conj(V(:,vterms));

%% OTHER FUNCTIONS

function v = myGetField(s,f), v = s; for i = 1:numel(f), v = v.(f{i}); end

function p = uifunction(w,h,X,t)
x = X(:,1);
y = X(:,2);
k = computeWaveNumber(w,h);
p = exp(1i*k*(x*cos(t) + y*sin(t)));



