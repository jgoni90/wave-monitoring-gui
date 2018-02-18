function S = PGDmapping_pg_K(idim,mode,matrices,coefs,parameters,rhs,algorithm,pgdRB)

%% IMPLEMENTATION OF THE PGD MAPPING S = S(R,T,...,Fn,Gn,Hn,...)

restdims = 1:parameters.nOfPGDdimensions;
restdims(idim) = [];
restnames = parameters.PGDdimensions(restdims,1);
im = sqrt(-1);

%Implementations
if mystrcmp('XY',restnames,true)
    
    %Some definitions
    R  = mode{restdims}(:,1);
    Rd = mode{restdims}(:,2);
    Rdt = Rd';
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
            iRt = Rdt * matrices(restdims).(ifields{IC}){i};
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
            alphaR = Rdt * rhsXY{j};
            iF = iF + rhsK{j} * alphaR.';
        end
        prevA = prevA + icoef*iF;    
    end
    
    %Compute mode S(w)
    [lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
    S = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));
    
    %Compute dual Ss
    iS = matrices(idim).M * S;
    fs = (R'*matrices(restdims).M*R)' * iS;
    Ss = lu_P'*(lu_L'\(lu_U'\(lu_Q'*fs)));
    
%------------
    
elseif mystrcmp({'XY' 'THETA'},restnames,true)
    
    %Dimensions
    restdimsaux   = findPGDdimension({'XY' 'THETA'},parameters.PGDdimensions);
    dimXY         = restdimsaux(1);
    dimTHETA      = restdimsaux(2);
    nOfnodesK     = size(matrices(idim).M,1);
    nOfnodesTHETA = size(matrices(dimTHETA).M,1);
    nOfrhsTerms   = size(rhs{dimXY}{1}{1},2);
    
    %Some definitions
    R  = mode{dimXY}(:,1);
    Rd = mode{dimXY}(:,2);
    Rdt = Rd';
    T  = mode{dimTHETA}(:,1);
    Td = mode{dimTHETA}(:,2);
    Tdt = Td';
    iTt = Tdt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    Tc_reshaped = reshape(Tdt,1,nOfnodesTHETA,1);
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
            iRt = Rdt * matrices(dimXY).(ifields{IC}){i};
            A = A + icoefs(IC)*(iRt*R) * matrices(idim).varM{coefIC}{i};
            iF = iRt * pgdRB{dimXY};
            prevA = prevA - matrices(idim).varM{coefIC}{i} * (pgdRB{idim}*(icoefs(IC)*iF.*iH).');
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;

    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        rhsKTHETA = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsKTHETA)
            alphaR = Rdt * rhsXY{j};
            ijF = bsxfun(@times,reshape(rhsKTHETA{j},nOfnodesK,nOfnodesTHETA,nOfrhsTerms),Tc_reshaped);
            ijF = reshape(sum(ijF,2),nOfnodesK,nOfrhsTerms);
            iF = iF + ijF * alphaR.';
        end
        prevA = prevA + icoef*iF;    
    end
    
    %Compute mode S(w)
    [lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
    S = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));
    
    %Compute dual Ss
    iS = matrices(idim).M * S;
    fs = ((R'*matrices(dimXY).M*R) * (T'*matrices(dimTHETA).M*T))' * iS;
    Ss = lu_P'*(lu_L'\(lu_U'\(lu_Q'*fs)));
    
%------------  %Hereafter are the implementations with ALPHA*. Use first the ones with the larger number of dimensions

elseif mystrcmp({'XY' 'THETA' 'ALPHA*'},restnames)
    
    %Dimensions
    restdimsaux   = findPGDdimension({'XY' 'THETA'},parameters.PGDdimensions);
    dimXY         = restdimsaux(1);
    dimTHETA      = restdimsaux(2);
    dimALPHA      = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA        = length(dimALPHA);
    nOfnodesK     = size(matrices(idim).M,1);
    nOfnodesTHETA = size(matrices(dimTHETA).M,1);
    nOfrhsTerms   = size(rhs{dimXY}{1}{1},2);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    R  = mode{dimXY}(:,1);
    Rd = mode{dimXY}(:,2);
    Rdt = Rd';
    T  = mode{dimTHETA}(:,1);
    Td = mode{dimTHETA}(:,2);
    Tdt = Td';
    iTt = Tdt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    Tc_reshaped = reshape(Tdt,1,nOfnodesTHETA,1);
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
    iALPHAdualcoef = 1;
    iP = ones(1,nOfprevTerms);
    iPtcoef = zeros(1,nALPHA); ivarPtcoef = iPtcoef;
    ViPtcoef = zeros(nALPHA,nOfprevTerms); VivarPtcoef = ViPtcoef;
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        P  = mode{idimALPHA}(:,1);
        Pd = mode{idimALPHA}(:,2);
        Pdt = Pd';
        iPt = Pdt * matrices(idimALPHA).M;
        ivarPt = Pdt * matrices(idimALPHA).varM;
        
        %Store
        iPtcoef(iALPHA)       = iPt * P;
        ViPtcoef(iALPHA,:)    = iPt * pgdRB{idimALPHA};
        ivarPtcoef(iALPHA)    = ivarPt * P;
        VivarPtcoef(iALPHA,:) = ivarPt * pgdRB{idimALPHA};
        
        %Coefficients for the weak form terms without ALPHA boundaries
        iALPHAcoef     = iALPHAcoef * iPtcoef(iALPHA);
        iP             = iP .* ViPtcoef(iALPHA,:);
        iALPHArhscoef  = iALPHArhscoef * (Pdt*rhs{idimALPHA}); %for RHS
        iALPHAdualcoef = iALPHAdualcoef * (P'*matrices(idimALPHA).M*P); %for dual
    end
    iPiH = iP.*iH;
    
    %Weak form terms without ALPHA boundaries
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rdt * matrices(dimXY).(ifields{IC}){i};
            A = A + icoefs(IC)*iALPHAcoef*(iRt*R) * matrices(idim).varM{coefIC}{i};
            iF = iRt * pgdRB{dimXY};
            prevA = prevA - matrices(idim).varM{coefIC}{i} * (pgdRB{idim}*(icoefs(IC)*iF.*iPiH).');
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
            iRt = Rdt * matrices(dimXY).C_alpha{idimALPHA}{i};
            A = A + im*ivarALPHAcoef(iALPHA)*(iRt*R) * matrices(idim).varM{coefsINDEX(end)}{i};
            iF = iRt * pgdRB{dimXY};
            prevA = prevA - matrices(idim).varM{coefsINDEX(end)}{i} * (pgdRB{idim}*(im*iF.*ivarP(iALPHA,:).*iH).');
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        rhsKTHETA = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsKTHETA)
            alphaR = Rdt * rhsXY{j};
            ijF = bsxfun(@times,reshape(rhsKTHETA{j},nOfnodesK,nOfnodesTHETA,nOfrhsTerms),Tc_reshaped);
            ijF = reshape(sum(ijF,2),nOfnodesK,nOfrhsTerms);
            iF = iF + ijF * alphaR.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
    
    %Compute mode S(w)
    [lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
    S = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));
    
    %Compute dual Ss
    iS = matrices(idim).M * S;
    fs = (iALPHAdualcoef * (R'*matrices(dimXY).M*R) * (T'*matrices(dimTHETA).M*T))' * iS;
    Ss = lu_P'*(lu_L'\(lu_U'\(lu_Q'*fs)));
    
%------------

elseif mystrcmp({'XY' 'ALPHA*'},restnames)
    
    %Dimensions
    dimXY     = findPGDdimension('XY',parameters.PGDdimensions);
    dimALPHA  = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA    = length(dimALPHA);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    R  = mode{dimXY}(:,1);
    Rd = mode{dimXY}(:,2);
    Rdt = Rd';
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
    iALPHAdualcoef = 1;
    iP = ones(1,nOfprevTerms);
    iPtcoef = zeros(1,nALPHA); ivarPtcoef = iPtcoef;
    ViPtcoef = zeros(nALPHA,nOfprevTerms); VivarPtcoef = ViPtcoef;
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        P  = mode{idimALPHA}(:,1);
        Pd = mode{idimALPHA}(:,2);
        Pdt = Pd';
        iPt = Pdt * matrices(idimALPHA).M;
        ivarPt = Pdt * matrices(idimALPHA).varM;
        
        %Store
        iPtcoef(iALPHA)       = iPt * P;
        ViPtcoef(iALPHA,:)    = iPt * pgdRB{idimALPHA};
        ivarPtcoef(iALPHA)    = ivarPt * P;
        VivarPtcoef(iALPHA,:) = ivarPt * pgdRB{idimALPHA};
        
        %Coefficients for the weak form terms without ALPHA boundaries
        iALPHAcoef     = iALPHAcoef * iPtcoef(iALPHA);
        iP             = iP .* ViPtcoef(iALPHA,:);
        iALPHArhscoef  = iALPHArhscoef * (Pdt*rhs{idimALPHA}); %for RHS
        iALPHAdualcoef = iALPHAdualcoef * (P'*matrices(idimALPHA).M*P); %for dual
    end
    
    %Weak form terms without ALPHA boundaries
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rdt * matrices(dimXY).(ifields{IC}){i};
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
            iRt = Rdt * matrices(dimXY).C_alpha{idimALPHA}{i};
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
            alphaR = Rdt * rhsXY{j};
            iF = iF + rhsK{j} * alphaR.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
    
    %Compute mode S(w)
    [lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
    S = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));
    
    %Compute dual Ss
    iS = matrices(idim).M * S;
    fs = (iALPHAdualcoef * (R'*matrices(dimXY).M*R))' * iS;
    Ss = lu_P'*(lu_L'\(lu_U'\(lu_Q'*fs)));
end

%Normalize if needed
if idim < parameters.nOfPGDdimensions
    S2 = S' * iS;
    S = S/sqrt(S2);
    S2 = Ss' * matrices(idim).M * Ss;
    Ss = Ss/sqrt(S2);
end

%Store
S = [S,Ss];

