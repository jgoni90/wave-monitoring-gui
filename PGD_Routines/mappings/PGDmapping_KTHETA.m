function S = PGDmapping_KTHETA(idim,mode,matrices,coefs,parameters,rhs,algorithm,pgdRB)

%% IMPLEMENTATION OF THE PGD MAPPING S = S(R,T,...,Fn,Gn,Hn,...)

restdims = 1:parameters.nOfPGDdimensions;
restdims(idim) = [];
restnames = parameters.PGDdimensions(restdims,1);
im = sqrt(-1);

%Implementations
if mystrcmp('XY',restnames,true)
    
    %Some definitions
    R = mode{restdims};
    Rt = R';
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    nOfWeakTerms = 6;
    ifields = {'Kint','Kxpml','Kypml','varM','C','D'};
    coefsINDEX = [1:5 5];
    icoefs  = [-1,-1,-1,1,im,im];
    
    %Weak form terms
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(restdims).(ifields{IC}){i};
            A = A + icoefs(IC)*(iRt*R) * matrices(idim).varM{coefIC}{i};
            iF = iRt * pgdRB{restdims};
            prevA = prevA - matrices(idim).varM{coefIC}{i} * (pgdRB{idim}*(icoefs(IC)*iF.'));
        end
    end
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{restdims}{i};
        rhsK = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsK)
            alphaR = Rt * rhsXY{j};
            iF = iF + rhsK{j} * alphaR.';
        end
        prevA = prevA + icoef*iF;    
    end
    
%------------ %Hereafter are the implementations with ALPHA*. Use first the ones with the larger number of dimensions

elseif mystrcmp({'XY' 'ALPHA*'},restnames)
    
    %Dimensions
    dimXY     = findPGDdimension('XY',parameters.PGDdimensions);
    dimALPHA  = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA    = length(dimALPHA);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    R = mode{dimXY};
    Rt = R';
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    if isempty(matrices(dimXY).C)
        nOfWeakTerms = 5;
        ifields = {'Kint','Kxpml','Kypml','varM','D'};
        coefsINDEX = 1:5;
        icoefs  = [-1,-1,-1,1,im];
    else
        nOfWeakTerms = 6;
        ifields = {'Kint','Kxpml','Kypml','varM','C','D'};
        coefsINDEX = [1:5 5];
        icoefs = [-1,-1,-1,1,im,im];
    end
    
    %Precomputations for the ALPHA dimensions for speed
    iALPHAcoef = 1;
    iALPHArhscoef = 1;
    iP = ones(1,nOfprevTerms);
    iPtcoef = zeros(1,nALPHA); ivarPtcoef = iPtcoef;
    ViPtcoef = zeros(nALPHA,nOfprevTerms); VivarPtcoef = ViPtcoef;
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        P = mode{idimALPHA};
        Pt = P';
        iPt = Pt * matrices(idimALPHA).M;
        ivarPt = Pt * matrices(idimALPHA).varM;
        
        %Store
        iPtcoef(iALPHA)       = iPt * P;
        ViPtcoef(iALPHA,:)    = iPt * pgdRB{idimALPHA};
        ivarPtcoef(iALPHA)    = ivarPt * P;
        VivarPtcoef(iALPHA,:) = ivarPt * pgdRB{idimALPHA};
        
        %Coefficients for the weak form terms without ALPHA boundaries
        iALPHAcoef    = iALPHAcoef * iPtcoef(iALPHA);
        iP            = iP .* ViPtcoef(iALPHA,:);
        iALPHArhscoef = iALPHArhscoef * (Pt*rhs{idimALPHA}); %for RHS
    end
    
    %Weak form terms without ALPHA boundaries
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            A = A + icoefs(IC)*iALPHAcoef*(iRt*R) * matrices(idim).varM{coefIC}{i};
            iF = iRt * pgdRB{dimXY};
            prevA = prevA - matrices(idim).varM{coefIC}{i} * (pgdRB{idim}*(icoefs(IC)*iF.*iP).');
        end
    end
    
    %Coefficients for the weak form terms with ALPHA boundaries
    ivarALPHAcoef = ivarPtcoef;
    ivarP = VivarPtcoef;
    for iALPHA = 1:nALPHA
        for jALPHA = setdiff(1:nALPHA,iALPHA)
            ivarALPHAcoef(iALPHA) = ivarALPHAcoef(iALPHA) * iPtcoef(jALPHA);
            ivarP(iALPHA,:) = ivarP(iALPHA,:) .* ViPtcoef(jALPHA,:);
        end
    end
    
    %Add ALPHA boundaries (using the last coefficient)
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            iRt = Rt * matrices(dimXY).C_alpha{idimALPHA}{i};
            A = A + im*ivarALPHAcoef(iALPHA)*(iRt*R) * matrices(idim).varM{coefsINDEX(end)}{i};
            iF = iRt * pgdRB{dimXY};
            prevA = prevA - matrices(idim).varM{coefsINDEX(end)}{i} * (pgdRB{idim}*(im*iF.*ivarP(iALPHA,:)).');
        end
    end
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        rhsK = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsK)
            alphaR = Rt * rhsXY{j};
            iF = iF + rhsK{j} * alphaR.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
end

%Compute mode S(w)
[lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
S = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));

%Normalize if needed
if idim < parameters.nOfPGDdimensions
    S2 = S' * matrices(idim).M * S;
    S = S/sqrt(S2);
end

