function coefs = computeSVDcoefs_XY_K(parameters,PGDmeshes)

%% DEFININITION OF COEFFICIENTS

dimXY = findPGDdimension('XY',parameters.PGDdimensions);
dimK  = findPGDdimension({'K' 'KTHETA'},parameters.PGDdimensions);
dimK  = dimK(logical(dimK));

%coef1 = ccg in \Omega_int x I_omega
coefs(1).handle         = @coef1;
coefs(1).extraInputs    = {};
coefs(1).PGDmeshesField = {'T','int'};

%coef2 = ccg * (sigma_y/sigma_x) in \Omega_pml x I_omega
coefs(2).handle         = @coef2;
coefs(2).extraInputs    = {{'PML','sigma'}};
coefs(2).PGDmeshesField = {'T','ext'};

%coef3 = ccg * (sigma_x/sigma_y) in \Omega_pml x I_omega
coefs(3).handle         = @coef3;
coefs(3).extraInputs    = {{'PML','sigma'}};
coefs(3).PGDmeshesField = {'T','ext'};

%coef4 = ccg * sigma_x * sigma_y * k^2 in \Omega x I_omega
coefs(4).handle         = @coef4;
coefs(4).extraInputs    = {{'PML','sigma'}};
coefs(4).PGDmeshesField = {'T','all'};

%coef5 = ccg * k in (\Gamma_R \cup \Gamma_pml) x I_omega
coefs(5).handle         = @coef5;
coefs(5).extraInputs    = {};
coefs(5).PGDmeshesField = {'Tb','all'};

%% COMPUTATION OF SVD

nOfCoefs = numel(coefs);
tolSVD = parameters.coefparam.relerror;
nOfSVDfunctions = 2;
nOfNodes_xy = size(PGDmeshes(dimXY).X,1);
onesxy = ones(nOfNodes_xy,1);
nOfNodes_w = size(PGDmeshes(dimK).X,1);
for i = 1:nOfCoefs
    
    disp(['  Coefficient ' num2str(i) ':'])
    
    %Spatial nodes
    T = myGetField(PGDmeshes(dimXY),coefs(i).PGDmeshesField);
    xy_nodes = unique(T);
    
    %Inputs for computing the coefficients
    nOfExtraInputs = numel(coefs(i).extraInputs);
    extraInputs = cell(nOfExtraInputs,1);
    for j = 1:nOfExtraInputs
        v = myGetField(parameters.meshes(dimXY),coefs(i).extraInputs{j});
        extraInputs{j} = v(xy_nodes,:);
    end
    
    %SVD data
    disp('      Computation of matrix data...')
    SVDdata = zeros(size(xy_nodes,1),nOfNodes_w);
    for j = 1:nOfNodes_w
        SVDdata(:,j) = coefs(i).handle(PGDmeshes(dimK).X(j),...
                                       parameters.meshes(dimXY).bottom(xy_nodes),...
                                       extraInputs);
    end

    %SVD calculation
    disp('      Computation of SVD...')
    [U,s,V] = svd(SVDdata,0); %economic size SVD
    sdiag = diag(s);
    nOfTerms = find(sdiag < tolSVD,1);
    vterms = 1:nOfTerms;
    R = U(:,vterms)*s(vterms,vterms)*V(:,vterms)';
    err = norm(SVDdata(:) - R(:)) / norm(R(:));

    %Storage (only necessary terms)
    disp(['      ' num2str(nOfTerms) ' terms for an alpha tolerance ' num2str(tolSVD) ' with error ' num2str(err)])
    coefs(i).nOfTerms = nOfTerms;
    coefs(i).alpha = sdiag(vterms).';
    coefs(i).RB = cell(nOfSVDfunctions,1);
    F = zeros(nOfNodes_xy,nOfTerms);
    F(xy_nodes,:) = U(:,vterms);
    coefs(i).RB{1} = coefs(i).alpha(onesxy,:) .* F;
    coefs(i).RB{2} = conj(V(:,vterms));
end

%% OTHER FUNCTIONS

function v = myGetField(s,f), v = s; for i = 1:numel(f), v = v.(f{i}); end

function p = coef1(w,h,varargin)
k = computeWaveNumber(w,h);
p = celerities(w,k,h);

function p = coef2(w,h,varargin)
sigma = varargin{1}{1};
k = computeWaveNumber(w,h);
p = celerities(w,k,h) .* ((w+1i*sigma(:,2))./(w+1i*sigma(:,1)));

function p = coef3(w,h,varargin)
sigma = varargin{1}{1};
k = computeWaveNumber(w,h);
p = celerities(w,k,h) .* ((w+1i*sigma(:,1))./(w+1i*sigma(:,2)));

function p = coef4(w,h,varargin)
sigma = varargin{1}{1};
k = computeWaveNumber(w,h);
p = celerities(w,k,h) .* (1+ 1i*sigma(:,1)/w) .* (1+ 1i*sigma(:,2)/w) .* (k.^2);

function p = coef5(w,h,varargin)
k = computeWaveNumber(w,h);
p = celerities(w,k,h) .* k;











