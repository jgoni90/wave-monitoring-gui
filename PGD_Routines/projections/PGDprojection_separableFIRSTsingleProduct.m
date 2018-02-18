function projectedpgd = PGDprojection_separableFIRSTsingleProduct...
    (pgd,matrices,meshes,rhs,parameters,algorithm,maxGBmemory)

%% IMPLEMENTATION OF L2 PROJECTION: find u, separable, s.t. ||upgd - u||^2 = 0

%This function uses an standard PGD algorithm for solving the projection.
%upgd.RB is defined as upgd = pgd{1} * pgd{2}.
%It is defined ONLY for two separated dimensions.
%The second upgd dimension pgd{1}.RB{2} * pgd{2}.RB{2} is precomputed.
%The first upgd dimension pgd{1}.RB{1} * pgd{2}.RB{1} is vectorized according to maxGBmemory limit.

%Initialize projectedpgd structure and parameters
algorithmcopy                   = algorithm;
parameterscopy                  = parameters;
algorithmcopy.dual.value        = false;
algorithmcopy.update.value      = false;
algorithmcopy.projection.value  = false;
parameterscopy.recoverPGD.value = false;
parameterscopy.residualType     = algorithm.projection.residualType;
projectedpgd = initializePGD(algorithmcopy,meshes,matrices,parameterscopy);

%Initialize general variables
Uconvergence = false;
mode = cell(parameters.nOfPGDdimensions,1);
coefs = zeros(parameters.nOfPGDdimensions,1);
nnodes = zeros(parameters.nOfPGDdimensions,1);
coefprev = 1;
for idim = 1:parameters.nOfPGDdimensions
    nnodes(idim) = size(projectedpgd.RB{idim},1);
    mode{idim} = ones(size(projectedpgd.RB{idim}(:,1)));
    coefs(idim) = mode{idim}' * matrices(idim).M * mode{idim};
    coefprev = coefprev * coefs(idim);
end
previousMode = mode;
prevresidual = 0;
if strcmpi(parameterscopy.residualType,'DRFS')
    TOL_proj_residual = (pgd.errors.residual{pgd.counters.projTerm}(pgd.counters.proj-1) +...
        pgd.errors.residual{pgd.counters.projTerm}(pgd.counters.proj))/2;
end
MatCoeff = 0;
fCoeff = 0;
projectedpgd.alpha = 0;

%Copy of the original pgd for speed
nOfProducts       = 1;
pgdRB             = cell(nOfProducts+1,1);
precomputedTensor = cell(nOfProducts+1,1);
for j = 1:nOfProducts+1
    POSf  = 1:pgd{j}.counters.Uterm;
    pgdRB{j} = cell(1,parameters.nOfPGDdimensions);
    for i = 1:parameters.nOfPGDdimensions
        pgdRB{j}{i} = pgd{j}.RB{i}(:,POSf);
    end
end

%Complete precomputed tensor for dimension 2
nterms1 = size(pgdRB{1}{1},2);
nterms2 = size(pgdRB{2}{1},2);
onesrv  = ones(1,nterms2);
precomputedTensor{2} = bsxfun(@times,reshape(pgdRB{1}{2},nnodes(2),nterms1,1),...
    reshape(pgdRB{2}{2},nnodes(2),1,nterms2));

%Partial precomputed tensor for dimension 1 (it depedens on the max memory allowed by user)
auxnumber     = pgdRB{1}{1}(1);
whosauxnumber = whos('auxnumber');
if whosauxnumber.complex, coefcomplex = 2; else coefcomplex = 1; end
restGBmemory  = maxGBmemory - MSL_getGb(precomputedTensor{2});
ntermss       = floor(restGBmemory*(1024^3) / (8*coefcomplex*nnodes(1)*nterms2));
if ntermss > nterms1, ntermss = nterms1; end
prestoredPOS  = 1:ntermss;
precomputedTensor{1} = bsxfun(@times,reshape(pgdRB{1}{1}(:,prestoredPOS),nnodes(1),ntermss,1),...
    reshape(pgdRB{2}{1},nnodes(1),1,nterms2));

dim1Name = parameters.PGDdimensions{1,1};
disp([num2str(ntermss) ' terms in dimension ' dim1Name ...
    ' for a max memory allowed of ' num2str(restGBmemory) ' Gb'])

while ~Uconvergence && projectedpgd.counters.Uterm < parameterscopy.maxterms
    
    PGDconvergence = false;
    iter = 0;
    POSfn = 1:projectedpgd.counters.Uterm;

    while ~PGDconvergence && iter < algorithm.projection.maxiter
        
        normcoef = 1;
        normcoefprev = 1;

        for idim = 1:parameters.nOfPGDdimensions

            %Get dimensions
            restdims = 1:parameters.nOfPGDdimensions;
            restdims(idim) = [];

            %Computation of vectors for dimension restdims
            vM = mode{restdims}'*matrices(restdims).M;

            %Coefficient for the system matrix A
            acoef = coefs(restdims);

            %Coefficient (vectorized) for the system vector fn
            vMfn = vM*projectedpgd.RB{restdims}(:,POSfn);
            
            %Product with alpha coefficients to the previous reduced basis
            vMfn = vMfn .* projectedpgd.alpha;
            
            %Vector f
            if idim == 1
                %Precomputed terms
                vMf = vM * reshape(precomputedTensor{restdims},nnodes(restdims),nterms1*nterms2);
                vMf = reshape(vMf,nterms1,nterms2);
                prevMf = reshape(vMf(prestoredPOS,:),ntermss*nterms2,1);
                fcoef  = reshape(precomputedTensor{idim},nnodes(idim),ntermss*nterms2) * prevMf;
                
                %Non-precomputed terms
                for iterm = ntermss+1:nterms1
                    iprevMf = reshape(vMf(iterm,:),nterms2,1);
                    itermmode = pgdRB{1}{idim}(:,iterm);
                    fcoef = fcoef + (itermmode(:,onesrv) .* pgdRB{2}{idim}) * iprevMf;
                end
            else
                %Precomputed terms
                prevMf = vM * reshape(precomputedTensor{restdims},nnodes(restdims),ntermss*nterms2);
                fcoef  = reshape(precomputedTensor{idim}(:,prestoredPOS,:),nnodes(idim),ntermss*nterms2) * prevMf.';
                
                %Non-precomputed terms
                for iterm = ntermss+1:nterms1
                    itermmode = pgdRB{1}{restdims}(:,iterm);
                    vMf = vM * (itermmode(:,onesrv) .* pgdRB{2}{restdims});
                    itermmode2 = reshape(precomputedTensor{idim}(:,iterm,:),nnodes(idim),nterms2);
                    fcoef = fcoef + itermmode2 * vMf.';
                end
            end
            
            %Vector fn
            if projectedpgd.counters.Uterm > 0
                fncoef = projectedpgd.RB{idim}(:,POSfn)*vMfn.';
            else
                fncoef = 0;
            end

            %Compute the mode for the current PGD dimension
            F = (1/acoef)*(fcoef - fncoef);
            
            %Normalize if needed
            if ~algorithm.projection.projCoeff && idim > 1
                modenorm = F' * matrices(idim).M * F; %normalize all dimensions but the first one
                F = F/sqrt(modenorm);
            elseif algorithm.projection.projCoeff 
                modenorm = F' * matrices(idim).M * F; %normalize all dimensions
                F = F/sqrt(modenorm);
            end

            %Some updates
            coefs(idim) = F' * matrices(idim).M * F;
            normcoef = normcoef * coefs(idim);
            normcoefprev = normcoefprev * F'*matrices(idim).M*previousMode{idim};
            mode{idim} = F;
        end

        %Error and convergence
        PGDabserror = sqrt(real(normcoef + coefprev - normcoefprev - normcoefprev'));
        PGDerror = PGDabserror/sqrt(real(normcoef));
        PGDconvergence = PGDerror < parameters.toliter;

        %Update
        previousMode = mode; 
        coefprev = normcoef;
        iter = iter + 1;
    end
    
    %Compute projection of coefficients if needed
    if algorithm.projection.projCoeff
        [projectedpgd.alpha,MatCoeff,fCoeff] = computeProjectionCoefficients...
            (projectedpgd,mode,MatCoeff,fCoeff,matrices,parameters,...
            pgdRB,precomputedTensor,nnodes,ntermss,nterms1,nterms2,prestoredPOS,onesrv);
    else
        projectedpgd.alpha = ones(1,projectedpgd.counters.Uterm+1);
    end
    
    %Update term
    projectedpgd.counters.Uterm = projectedpgd.counters.Uterm + 1;
    for idim = 1:parameters.nOfPGDdimensions
        projectedpgd.RB{idim}(:,projectedpgd.counters.Uterm) = mode{idim};
    end
    
    %Compute error in U
    if strcmpi(parameterscopy.residualType,'ALPHACOEF')
        residual = norm(projectedpgd.alpha(end)) / norm(projectedpgd.alpha);
    else
        recursiveComputation = ~algorithm.projection.projCoeff;
        [residual,projectedpgd.errors.residualMat] = ...
                computePGDerror(parameterscopy,projectedpgd,meshes,matrices,rhs,algorithm,recursiveComputation);
    end
    projectedpgd.counters.residualUterm = projectedpgd.counters.residualUterm + 1;
    projectedpgd.errors.residual{1}(projectedpgd.counters.Uterm) = residual;
    
    %Convergence
    if strcmpi(parameterscopy.residualType,'DRFS')
        relerror = abs((residual - prevresidual)/residual);
        Uconvergence = relerror < algorithm.projection.relerror && ...
                       residual < TOL_proj_residual*algorithm.projection.tolfactor;
        prevresidual = residual;
    else
        Uconvergence = residual < algorithm.projection.relerror;
    end
    
    %Print
    printPGD(projectedpgd.counters.Uterm,residual,iter,PGDerror);

end

%Prepare output
projectedpgd.errors.residual{1}(projectedpgd.counters.Uterm+1:end) = [];
for idim = 1:parameters.nOfPGDdimensions
    projectedpgd.RB{idim}(:,projectedpgd.counters.Uterm+1:end) = [];
end

%Add alpha coefficients to the final reduced basis
if algorithm.projection.projCoeff
    auxones = ones(size(matrices(1).M,1),1);
    projectedpgd.RB{1} = projectedpgd.RB{1} .* projectedpgd.alpha(auxones,:);
end

%% Projection function of coefficients

function [alpha,Mat,f] = computeProjectionCoefficients...
            (pgd,mode,prevMat,prevf,matrices,parameters,...
            pgdRB,precomputedTensor,nnodes,ntermss,nterms1,nterms2,prestoredPOS,onesrv)

%Initialize
Uterm                   = pgd.counters.Uterm;        
Mat                     = ones(Uterm+1,Uterm+1); %diag(Mat) = (mode,mode)_L2 = 1
f                       = ones(Uterm+1,1);
f(1:Uterm)              = prevf;
Mat(1:Uterm,1:Uterm)    = prevMat;
vstored                 = cell(1,parameters.nOfPGDdimensions);
for j = 1:parameters.nOfPGDdimensions, vstored{j} = matrices(j).M * mode{j}; end

%Compute last column -but diagonal term- of matrix Mat as (mode,pgd.RB(:)(:,1:Uterm))_L2
for i = 1:Uterm
    for j = 1:parameters.nOfPGDdimensions
        Mat(i,Uterm+1) = Mat(i,Uterm+1) * (pgd.RB{j}(:,i)' * vstored{j});
    end
end

%Complete the hermitian system matrix Mat
Mat(Uterm+1,1:Uterm) = Mat(1:Uterm,Uterm+1)';

%Last term of separable RHS vector computed with pgdRB structure
vM = cell(parameters.nOfPGDdimensions,1);
for j = 1:parameters.nOfPGDdimensions, vM{j} = mode{j}'*matrices(j).M; end

%Precomputed terms
v11 = vM{1} * reshape(precomputedTensor{1},nnodes(1),ntermss*nterms2);
v12 = vM{2} * reshape(precomputedTensor{2}(:,prestoredPOS,:),nnodes(2),ntermss*nterms2);
v1 = sum(v11.*v12);

%Non-precomputed terms
v2 = 0;
for iterm = ntermss+1:nterms1
    itermmode1 = pgdRB{1}{1}(:,iterm);
    v21 = vM{1} * (itermmode1(:,onesrv) .* pgdRB{2}{1});
    itermmode2 = reshape(precomputedTensor{2}(:,iterm,:),nnodes(2),nterms2);
    v22 = vM{2} * itermmode2;
    v2 = v2 + v21.*v22;
end
f(Uterm+1) = v1 + sum(v2);

%Compute alpha coefficients
alpha = (Mat \ f).';

