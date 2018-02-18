function projectedpgd = PGDprojection_separableSingleProduct(pgd,matrices,meshes,rhs,parameters,algorithm)

%% IMPLEMENTATION OF L2 PROJECTION: find u, separable, s.t. ||upgd - u||^2 = 0

%This function uses an standard PGD algorithm for solving the projection.
%upgd is defined as upgd = pgd{1} * pgd{2}.
%It is defined for any number of separated dimensions.

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
coefprev = 1;
vM = cell(parameters.nOfPGDdimensions,1);
for idim = 1:parameters.nOfPGDdimensions
    mode{idim} = ones(size(projectedpgd.RB{idim}(:,1)));
    coefs(idim) = mode{idim}' * matrices(idim).M * mode{idim};
    coefprev = coefprev * coefs(idim);
    vM{idim} = zeros(1,size(meshes(idim).X,1));
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
nOfProducts = 1;
POSf        = cell(nOfProducts+1,1);
pgdRB       = cell(nOfProducts+1,1);
for j = 1:nOfProducts+1
    POSf{j}  = 1:pgd{j}.counters.Uterm;
    pgdRB{j} = cell(1,parameters.nOfPGDdimensions);
    for i = 1:parameters.nOfPGDdimensions
        pgdRB{j}{i} = pgd{j}.RB{i}(:,POSf{j});
    end
end
onespgd2 = ones(1,pgd{2}.counters.Uterm);

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

            %Computations for all dimensions but idim
            acoef = 1;
            vMfn = 1;
            for irest = restdims

                vM{irest} = mode{irest}'*matrices(irest).M;

                %Coefficient for the system matrix A
                acoef = acoef * coefs(irest);

                %Coefficient (vectorized) for the system vector fn
                vMfn = vMfn .* (vM{irest}*projectedpgd.RB{irest}(:,POSfn));
            end
            
            %Product with alpha coefficients to the previous reduced basis
            vMfn = vMfn .* projectedpgd.alpha;
            
            %Vector f (vectorized for the terms of pgd{2})
            fcoef = 0;
            for iterm = POSf{1}
                vMf = 1;
                for irest = restdims
                    irestmode = pgdRB{1}{irest}(:,iterm);
                    vMf = vMf .* (vM{irest}*(irestmode(:,onespgd2) .* pgdRB{2}{irest}));
                end
                itermmode = pgdRB{1}{idim}(:,iterm);
                fcoef = fcoef + (itermmode(:,onespgd2) .* pgdRB{2}{idim}) * vMf.';
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
            (projectedpgd,mode,MatCoeff,fCoeff,matrices,parameters,pgdRB,POSf,onespgd2);
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
            (pgd,mode,prevMat,prevf,matrices,parameters,pgdRB,POSf,onespgd2)

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
v = 0;
vM = cell(parameters.nOfPGDdimensions,1);
for j = 1:parameters.nOfPGDdimensions, vM{j} = mode{j}'*matrices(j).M; end
for iterm = POSf{1}
    vMf = 1;
    for j = 1:parameters.nOfPGDdimensions
        jmode = pgdRB{1}{j}(:,iterm);
        vMf = vMf .* (vM{j}*(jmode(:,onespgd2) .* pgdRB{2}{j}));
    end
    v = v + vMf;
end
f(Uterm+1) = sum(v);

%Compute alpha coefficients
alpha = (Mat \ f).';

