function R = PGDmapping_XY(idim,mode,matrices,coefs,parameters,rhs,algorithm,pgdRB)

%% IMPLEMENTATION OF THE PGD MAPPING R = R(S,T,...,Fn,Gn,Hn,...)

paramdims = findPGDdimension('PARAMETRICDIM',parameters.PGDdimensions);
paramnames = parameters.PGDdimensions(paramdims,1);
im = sqrt(-1);

%Implementations for matrix A and vector prevA
if mystrcmp('KTHETA',paramnames,true) || mystrcmp('K',paramnames,true)
    
    %Some definitions
    S = mode{paramdims};
    St = S';
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
            iSt = St * matrices(paramdims).varM{coefIC}{i};
            A = A + icoefs(IC)*(iSt*S) * matrices(idim).(ifields{IC}){i};
            iG = iSt * pgdRB{paramdims};
            prevA = prevA - matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iG.'));
        end
    end
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsK = rhs{paramdims}{i};
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            alphaS = St * rhsK{j};
            iF = iF + rhsXY{j} * alphaS.';
        end
        prevA = prevA + icoef*iF;    
    end

%------------
    
elseif mystrcmp('THETA',paramnames,true)
    
    %Some definitions
    T = mode{paramdims};
    Tt = T';
    iTt = Tt * matrices(paramdims).M;
    iH = iTt * pgdRB{paramdims};
    alphaT = Tt * rhs{paramdims};
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
            A = A + icoefs(IC) * matrices(idim).(ifields{IC}){i};
            prevA = prevA - matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iH.'));
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    iF = 0;
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        jF = 0;
        for j = 1:numel(rhsXY)
            jF = jF + rhsXY{j};
        end
        iF = iF + icoef*jF;    
    end
    prevA = prevA + iF * alphaT.';
    
%------------    
    
elseif mystrcmp({'K' 'THETA'},paramnames,true)
    
    %Dimensions
    restdimsaux = findPGDdimension({'K' 'THETA'},parameters.PGDdimensions);
    dimK        = restdimsaux(1);
    dimTHETA    = restdimsaux(2);
    
    %Some definitions
    S = mode{dimK};
    St = S';
    T = mode{dimTHETA};
    Tt = T';
    iTt = Tt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    STt = conj(S) * Tt;
    STt = (STt(:)).';
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
            iSt = St * matrices(dimK).varM{coefIC}{i};
            A = A + icoefs(IC)*(iSt*S) * matrices(idim).(ifields{IC}){i};
            iG = iSt * pgdRB{dimK};
            prevA = prevA - matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iG.*iH).');
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;
    
    %Rest of the RHS vector (K x THETA dimension for the parameters)
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsKTHETA = rhs{dimK}{i};
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            alphaST = STt * rhsKTHETA{j};
            iF = iF + rhsXY{j} * alphaST.';
        end
        prevA = prevA + icoef * iF;    
    end
    
%------------ Hereafter are the implementations with ALPHA*. Use first the ones with the larger number of dimensions

elseif mystrcmp({'K' 'THETA' 'ALPHA*'},paramnames)
    
    %Dimensions
    restdimsaux = findPGDdimension({'K' 'THETA'},parameters.PGDdimensions);
    dimK        = restdimsaux(1);
    dimTHETA    = restdimsaux(2);
    dimALPHA    = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA      = length(dimALPHA);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    S = mode{dimK};
    St = S';
    T = mode{dimTHETA};
    Tt = T';
    iTt = Tt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    STt = conj(S) * Tt;
    STt = (STt(:)).';
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    if isempty(matrices(idim).C)
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
    iPiH = iP.*iH;
    
    %Weak form terms without ALPHA boundaries
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        iStS = zeros(coefs(coefIC).nOfTerms,1); %For speed (only values for coefsINDEX(end) are needed)
        iGv = zeros(coefs(coefIC).nOfTerms,nOfprevTerms);
        for i = 1:coefs(coefIC).nOfTerms
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iStS(i) = iSt*S;
            A = A + icoefs(IC)*iALPHAcoef*iStS(i) * matrices(idim).(ifields{IC}){i};
            iGv(i,:) = iSt * pgdRB{dimK};
            prevA = prevA - (matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iGv(i,:).*iPiH).'));
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
            A = A + im*ivarALPHAcoef(iALPHA)*iStS(i) * matrices(idim).C_alpha{idimALPHA}{i};
            prevA = prevA - (matrices(idim).C_alpha{idimALPHA}{i} * (pgdRB{idim}*(im*iGv(i,:).*ivarP(iALPHA,:).*iH).'));
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;
    
    %Rest of the RHS vector (K x THETA dimension for the parameters)
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsKTHETA = rhs{dimK}{i};
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            alphaST = STt * rhsKTHETA{j};
            iF = iF + rhsXY{j} * alphaST.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
    
%------------

elseif mystrcmp({'K' 'ALPHA*'},paramnames) | mystrcmp({'KTHETA' 'ALPHA*'},paramnames) %#ok<OR2>
    
    %Dimensions
    dimK      = findPGDdimension({'K' 'KTHETA'},parameters.PGDdimensions);
    dimK      = dimK(logical(dimK));
    dimALPHA  = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA    = length(dimALPHA);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    S = mode{dimK};
    St = S';
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    if isempty(matrices(idim).C)
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
        iStS = zeros(coefs(coefIC).nOfTerms,1); %For speed (only values for coefsINDEX(end) are needed)
        iGv = zeros(coefs(coefIC).nOfTerms,nOfprevTerms);
        for i = 1:coefs(coefIC).nOfTerms
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iStS(i) = iSt*S;
            A = A + icoefs(IC)*iALPHAcoef*iStS(i) * matrices(idim).(ifields{IC}){i};
            iGv(i,:) = iSt * pgdRB{dimK};
            prevA = prevA - (matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iGv(i,:).*iP).'));
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
            A = A + im*ivarALPHAcoef(iALPHA)*iStS(i) * matrices(idim).C_alpha{idimALPHA}{i};
            prevA = prevA - (matrices(idim).C_alpha{idimALPHA}{i} * (pgdRB{idim}*(im*iGv(i,:).*ivarP(iALPHA,:)).'));
        end
    end
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsK = rhs{dimK}{i};
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            alphaS = St * rhsK{j};
            iF = iF + rhsXY{j} * alphaS.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
    
%------------

elseif mystrcmp({'THETA' 'ALPHA*'},paramnames)
    
    %Dimensions
    dimTHETA  = findPGDdimension('THETA',parameters.PGDdimensions);
    dimALPHA  = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA    = length(dimALPHA);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    T = mode{dimTHETA};
    Tt = T';
    iTt = Tt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    alphaT = Tt * rhs{dimTHETA};
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    if isempty(matrices(idim).C)
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
    iPiH = iP.*iH;
    
    %Weak form terms without ALPHA boundaries
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            A = A + icoefs(IC)*iALPHAcoef * matrices(idim).(ifields{IC}){i};
            prevA = prevA - (matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC).*iPiH).'));
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
            A = A + im*ivarALPHAcoef(iALPHA) * matrices(idim).C_alpha{idimALPHA}{i};
            prevA = prevA - (matrices(idim).C_alpha{idimALPHA}{i} * (pgdRB{idim}*(im*ivarP(iALPHA,:).*iH).'));
        end
    end
    
    %Add THETA coefficient (LHS of the problem does not depend on THETA)
    A = (iTt*T) * A;
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    iF = 0;
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        jF = 0;
        for j = 1:numel(rhsXY)
            jF = jF + rhsXY{j};
        end
        iF = iF + icoef*jF;    
    end
    prevA = prevA + iALPHArhscoef * (iF*alphaT.');   
    
%------------

elseif mystrcmp('ALPHA*',paramnames)
    
    %Dimensions
    nALPHA = length(paramdims);
    
    %Some definitions
    nOfprevTerms = size(pgdRB{1},2);
    nidim = size(matrices(idim).M,1);
    NNZ = nnz(matrices(idim).M);
    A = spalloc(nidim,nidim,NNZ);
    prevA = 0;
    if isempty(matrices(idim).C)
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
        idimALPHA = paramdims(iALPHA);
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
            A = A + icoefs(IC)*iALPHAcoef * matrices(idim).(ifields{IC}){i};
            prevA = prevA - (matrices(idim).(ifields{IC}){i} * (pgdRB{idim}*(icoefs(IC)*iP).'));
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
        idimALPHA = paramdims(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            A = A + im*ivarALPHAcoef(iALPHA) * matrices(idim).C_alpha{idimALPHA}{i};
            prevA = prevA - (matrices(idim).C_alpha{idimALPHA}{i} * (pgdRB{idim}*(im*ivarP(iALPHA,:)).'));
        end
    end
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{idim}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            iF = iF + rhsXY{j};
        end
        prevA = prevA + iALPHArhscoef*icoef*sum(iF,2);    
    end
end

%Compute mode R(x)
[lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
R = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));

%Normalize if needed
if idim < parameters.nOfPGDdimensions
    R2 = R' * matrices(idim).M * R;
    R = R/sqrt(R2);
end



