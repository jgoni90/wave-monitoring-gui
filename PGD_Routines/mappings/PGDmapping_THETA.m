function T = PGDmapping_THETA(idim,mode,matrices,coefs,parameters,rhs,algorithm,pgdRB)
                             
%% IMPLEMENTATION OF THE PGD MAPPING T = T(R,S,...,Fn,Gn,Hn,...)

restdims = 1:parameters.nOfPGDdimensions;
restdims(idim) = [];
restnames = parameters.PGDdimensions(restdims,1);
im = sqrt(-1);

%Implementations for matrix A and vector prevA
if mystrcmp('XY',restnames,true)
    
    %Some definitions and initializations
    R = mode{restdims};
    Rt = R';
    coefXY = 0;
    prevCoefXY = 0;
    nOfWeakTerms = 6;
    ifields = {'Kint','Kxpml','Kypml','varM','C','D'};
    coefsINDEX = [1:5 5];
    icoefs  = [-1,-1,-1,1,im,im];
    
    %Weak form terms
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(restdims).(ifields{IC}){i};
            iF  = iRt * pgdRB{restdims};
            coefXY = coefXY + icoefs(IC) * (iRt*R);
            prevCoefXY = prevCoefXY - icoefs(IC) * iF;
        end
    end
    
    %System matrix and part of the RHS vector due to previous terms
    A     = coefXY * matrices(idim).M;
    prevA = matrices(idim).M * (pgdRB{idim} * prevCoefXY.');
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    iF = 0;
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{restdims}{i};
        icoef = rhscoefs(i);
        jF = 0;
        for j = 1:numel(rhsXY)
            jF = jF + rhsXY{j};
        end
        iF = iF + icoef*jF;    
    end
    prevA = prevA + rhs{idim} * (Rt*iF).';
    
%------------
    
elseif mystrcmp({'XY' 'K'},restnames,true)
    
    %Dimensions
    restdimsaux   = findPGDdimension({'XY' 'K'},parameters.PGDdimensions);
    dimXY         = restdimsaux(1);
    dimK          = restdimsaux(2);
    nOfnodesK     = size(matrices(dimK).M,1);
    nOfnodesTHETA = size(matrices(idim).M,1);
    nOfrhsTerms   = size(rhs{dimXY}{1}{1},2);
    nreshaped     = nOfnodesTHETA * nOfrhsTerms;
    
    %Some definitions and initializations
    R = mode{dimXY};
    Rt = R';
    S = mode{dimK};
    St = S';
    coefXYK = 0;
    prevCoefXYK = 0;
    nOfWeakTerms = 6;
    ifields = {'Kint','Kxpml','Kypml','varM','C','D'};
    coefsINDEX = [1:5 5];
    icoefs  = [-1,-1,-1,1,im,im];
    
    %Weak form terms
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iF  = iRt * pgdRB{dimXY};
            iG  = iSt * pgdRB{dimK};
            coefXYK = coefXYK + icoefs(IC) * (iRt*R) * (iSt*S);
            prevCoefXYK = prevCoefXYK - icoefs(IC) * (iF.*iG);
        end
    end
    
    %System matrix and part of the RHS vector due to previous terms
    A     = coefXYK * matrices(idim).M;
    prevA = matrices(idim).M * (pgdRB{idim} * prevCoefXYK.');
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        rhsKTHETA = rhs{dimK}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsKTHETA)
            alphaR = Rt * rhsXY{j};
            ijF = St * reshape(rhsKTHETA{j},nOfnodesK,nreshaped);
            ijF = reshape(ijF,nOfnodesTHETA,nOfrhsTerms);
            iF = iF + ijF * alphaR.';
        end
        prevA = prevA + icoef*iF;    
    end
    
%------------ Hereafter are the implementations with ALPHA*. Use first the ones with the larger number of dimensions

elseif mystrcmp({'XY' 'K' 'ALPHA*'},restnames)
    
    %Dimensions
    restdimsaux   = findPGDdimension({'XY' 'K'},parameters.PGDdimensions);
    dimXY         = restdimsaux(1);
    dimK          = restdimsaux(2);
    dimALPHA      = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA        = length(dimALPHA);
    nOfnodesK     = size(matrices(dimK).M,1);
    nOfnodesTHETA = size(matrices(idim).M,1);
    nOfrhsTerms   = size(rhs{dimXY}{1}{1},2);
    nreshaped     = nOfnodesTHETA * nOfrhsTerms;
    
    %Some definitions and initializations
    nOfprevTerms = size(pgdRB{1},2);
    R = mode{dimXY};
    Rt = R';
    S = mode{dimK};
    St = S';
    coefXYK = 0;
    prevCoefXYKALPHA = 0;
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

    %Weak form terms
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        iStS = zeros(coefs(coefIC).nOfTerms,1); %For speed (only values for coefsINDEX(end) are needed)
        iGv = zeros(coefs(coefIC).nOfTerms,nOfprevTerms);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iStS(i) = iSt * S;
            iF  = iRt * pgdRB{dimXY};
            iGv(i,:) = iSt * pgdRB{dimK};
            coefXYK = coefXYK + icoefs(IC) * (iRt*R) * iStS(i);
            prevCoefXYKALPHA = prevCoefXYKALPHA - icoefs(IC) * (iF.*iGv(i,:).*iP);
        end
    end
    coefXYKALPHA = iALPHAcoef*coefXYK;
    
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
            iF  = iRt * pgdRB{dimXY};
            coefXYKALPHA = coefXYKALPHA + im*ivarALPHAcoef(iALPHA) * (iRt*R) * iStS(i);
            prevCoefXYKALPHA = prevCoefXYKALPHA - im * (iF.*iGv(i,:).*ivarP(iALPHA,:));
        end
    end
    
    %System matrix and part of the RHS vector due to previous terms
    A     = coefXYKALPHA * matrices(idim).M;
    prevA = matrices(idim).M * (pgdRB{idim} * prevCoefXYKALPHA.');
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        rhsKTHETA = rhs{dimK}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsKTHETA)
            alphaR = Rt * rhsXY{j};
            ijF = St * reshape(rhsKTHETA{j},nOfnodesK,nreshaped);
            ijF = reshape(ijF,nOfnodesTHETA,nOfrhsTerms);
            iF = iF + ijF * alphaR.';
        end
        prevA = prevA + iALPHArhscoef*icoef*iF;    
    end
    
%------------
      
elseif mystrcmp({'XY' 'ALPHA*'},restnames)
    
    %Dimensions
    dimXY    = findPGDdimension('XY',parameters.PGDdimensions);
    dimALPHA = findPGDdimension('ALPHA*',parameters.PGDdimensions);
    nALPHA   = length(dimALPHA);
    
    %Some definitions and initializations
    nOfprevTerms = size(pgdRB{1},2);
    R = mode{dimXY};
    Rt = R';
    coefXY = 0;
    prevCoefXYALPHA = 0;
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

    %Weak form terms
    for IC = 1:nOfWeakTerms
        coefIC = coefsINDEX(IC);
        for i = 1:coefs(coefIC).nOfTerms
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iF  = iRt * pgdRB{dimXY};
            coefXY = coefXY + icoefs(IC) * (iRt*R);
            prevCoefXYALPHA = prevCoefXYALPHA - icoefs(IC) * (iF.*iP);
        end
    end
    coefXYALPHA = iALPHAcoef * coefXY;
    
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
            iF  = iRt * pgdRB{dimXY};
            coefXYALPHA = coefXYALPHA + im*ivarALPHAcoef(iALPHA) * (iRt*R);
            prevCoefXYALPHA = prevCoefXYALPHA - im * (iF.*ivarP(iALPHA,:));
        end
    end
    
    %System matrix and part of the RHS vector due to previous terms
    A     = coefXYALPHA * matrices(idim).M;
    prevA = matrices(idim).M * (pgdRB{idim} * prevCoefXYALPHA.');
    
    %Rest of the RHS vector
    iF = 0;
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        icoef = rhscoefs(i);
        jF = 0;
        for j = 1:numel(rhsXY)
            jF = jF + rhsXY{j};
        end
        iF = iF + icoef*jF;   
    end
    prevA = prevA + iALPHArhscoef * (rhs{idim}*(Rt*iF).');
end

%Compute mode T(theta)
[lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
T = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));

%Normalize if needed
if idim < parameters.nOfPGDdimensions
    T2 = T' * matrices(idim).M * T;
    T = T/sqrt(T2);
end

