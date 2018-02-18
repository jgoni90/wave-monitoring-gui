function P = PGDmapping_ALPHA(idim,mode,matrices,coefs,parameters,rhs,algorithm,pgdRB)

%% IMPLEMENTATION OF THE PGD MAPPING P = P(R,S,T,...,Fn,Gn,Hn,...)

restdims  = 1:parameters.nOfPGDdimensions;
restdims(idim) = [];
restnames = parameters.PGDdimensions(restdims,1);
dimALPHA  = findPGDdimension('ALPHA*',parameters.PGDdimensions);
dimALPHA  = setdiff(dimALPHA,idim);
nALPHA    = length(dimALPHA);
im = sqrt(-1);

%Implementations for matrix A and vector prevA
%All the implementations include ALPHA*. Use first the ones with the larger number of dimensions
if mystrcmp({'XY' 'K' 'THETA'},restnames)
    
    %Dimensions
    restdimsaux = findPGDdimension({'XY' 'K' 'THETA'},parameters.PGDdimensions);
    dimXY       = restdimsaux(1);
    dimK        = restdimsaux(2);
    dimTHETA    = restdimsaux(3);
    
    %Some definitions and initializations
    R = mode{dimXY};
    Rt = R';
    S = mode{dimK};
    St = S';
    T = mode{dimTHETA};
    Tt = T';
    iTt = Tt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    STt = conj(S) * Tt;
    STt = (STt(:)).';
    nOfprevTerms = size(pgdRB{1},2);
    coefXYK_1 = 0; coefXYKALPHA_2 = 0; coefXYK_3 = 0;
    prevCoefXYKTALPHA_1 = 0; prevCoefXYKTALPHA_2 = 0; prevCoefXYKTALPHA_3 = 0;
    prevCoefXYKTALPHA_4 = 0;
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
        P      = mode{idimALPHA};
        Pt     = P';
        iPt    = Pt * matrices(idimALPHA).M;
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
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iStS(i) = iSt*S;
            iF = iRt * pgdRB{dimXY};
            iGv(i,:) = iSt * pgdRB{dimK};
            coefXYK_1 = coefXYK_1 + icoefs(IC) * (iRt*R) * iStS(i);
            prevCoefXYKTALPHA_1 = prevCoefXYKTALPHA_1 - icoefs(IC) * (iF.*iGv(i,:).*iPiH);
        end
    end
    coefXYKALPHA_1 = coefXYK_1 * iALPHAcoef;
    
    %Coefficients for the weak form terms with ALPHA boundaries
    ivarALPHAcoef = ivarPtcoef;
    ivarP = VivarPtcoef;
    for iALPHA = 1:nALPHA
        for jALPHA = setdiff(1:nALPHA,iALPHA)
            ivarALPHAcoef(iALPHA) = ivarALPHAcoef(iALPHA) * iPtcoef(jALPHA);
            ivarP(iALPHA,:) = ivarP(iALPHA,:) .* ViPtcoef(jALPHA,:);
        end
    end
    
    %Add ALPHA boundaries which are not idim (using the last coefficient)
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            iRt = Rt * matrices(dimXY).C_alpha{idimALPHA}{i};
            iF  = iRt * pgdRB{dimXY};
            coefXYKALPHA_2 = coefXYKALPHA_2 + im * ivarALPHAcoef(iALPHA) * (iRt*R) * iStS(i);
            prevCoefXYKTALPHA_2 = prevCoefXYKTALPHA_2 - im * (iF.*iGv(i,:).*ivarP(iALPHA,:).*iH);
        end
    end
    
    %Coefficient for ALPHA boundary corresponding to idim (using the last coefficient)
    for i = 1:coefs(coefsINDEX(end)).nOfTerms
        iRt = Rt * matrices(dimXY).C_alpha{idim}{i};
        iF  = iRt * pgdRB{dimXY};
        coefXYK_3 = coefXYK_3 + im * (iRt*R) * iStS(i);
        prevCoefXYKTALPHA_3 = prevCoefXYKTALPHA_3 - im * (iF.*iGv(i,:).*iPiH);
    end
    coefXYKALPHA_3 = coefXYK_3 * iALPHAcoef;
    
    %System matrix and part of the RHS vector due to previous terms
    coefTHETA = iTt * T;
    A     = coefTHETA*(coefXYKALPHA_1 + coefXYKALPHA_2)*matrices(idim).M + ...
            coefTHETA*coefXYKALPHA_3*matrices(idim).varM;
    prevA = matrices(idim).M * (pgdRB{idim} * (prevCoefXYKTALPHA_1 + prevCoefXYKTALPHA_2).') + ...
            matrices(idim).varM * (pgdRB{idim} * prevCoefXYKTALPHA_3.');
    
    %Rest of the RHS vector (K x THETA dimension for the parameters)
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsKTHETA = rhs{dimK}{i};
        rhsXY = rhs{dimXY}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            iF = iF + (Rt*rhsXY{j}) * (STt*rhsKTHETA{j}).';
        end
        prevCoefXYKTALPHA_4 = prevCoefXYKTALPHA_4 + icoef*iF;    
    end
    prevA = prevA + iALPHArhscoef*prevCoefXYKTALPHA_4*rhs{idim};
    
%------------    
    
elseif mystrcmp({'XY' 'K'},restnames) | mystrcmp({'XY' 'KTHETA'},restnames) %#ok<OR2>
    
    %Dimensions
    restdimsaux = findPGDdimension({'XY' 'K' 'KTHETA'},parameters.PGDdimensions);
    dimXY       = restdimsaux(1);
    dimK        = restdimsaux(2:3);
    dimK        = dimK(logical(dimK));
    
    %Some definitions and initializations
    R = mode{dimXY};
    Rt = R';
    S = mode{dimK};
    St = S';
    nOfprevTerms = size(pgdRB{1},2);
    coefXYK_1 = 0; coefXYKALPHA_2 = 0; coefXYK_3 = 0;
    prevCoefXYKALPHA_1 = 0; prevCoefXYKALPHA_2 = 0; prevCoefXYKALPHA_3 = 0;
    prevCoefXYKALPHA_4 = 0;
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
        P      = mode{idimALPHA};
        Pt     = P';
        iPt    = Pt * matrices(idimALPHA).M;
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
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iSt = St * matrices(dimK).varM{coefIC}{i};
            iStS(i) = iSt*S;
            iF = iRt * pgdRB{dimXY};
            iGv(i,:) = iSt * pgdRB{dimK};
            coefXYK_1 = coefXYK_1 + icoefs(IC) * (iRt*R) * iStS(i);
            prevCoefXYKALPHA_1 = prevCoefXYKALPHA_1 - icoefs(IC) * (iF.*iGv(i,:).*iP);
        end
    end
    coefXYKALPHA_1 = coefXYK_1 * iALPHAcoef;
    
    %Coefficients for the weak form terms with ALPHA boundaries
    ivarALPHAcoef = ivarPtcoef;
    ivarP = VivarPtcoef;
    for iALPHA = 1:nALPHA
        for jALPHA = setdiff(1:nALPHA,iALPHA)
            ivarALPHAcoef(iALPHA) = ivarALPHAcoef(iALPHA) * iPtcoef(jALPHA);
            ivarP(iALPHA,:) = ivarP(iALPHA,:) .* ViPtcoef(jALPHA,:);
        end
    end
    
    %Add ALPHA boundaries which are not idim (using the last coefficient)
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            iRt = Rt * matrices(dimXY).C_alpha{idimALPHA}{i};
            iF  = iRt * pgdRB{dimXY};
            coefXYKALPHA_2 = coefXYKALPHA_2 + im * ivarALPHAcoef(iALPHA) * (iRt*R) * iStS(i);
            prevCoefXYKALPHA_2 = prevCoefXYKALPHA_2 - im * (iF.*iGv(i,:).*ivarP(iALPHA,:));
        end
    end
    
    %Coefficient for ALPHA boundary corresponding to idim (using the last coefficient)
    for i = 1:coefs(coefsINDEX(end)).nOfTerms
        iRt = Rt * matrices(dimXY).C_alpha{idim}{i};
        iF  = iRt * pgdRB{dimXY};
        coefXYK_3 = coefXYK_3 + im * (iRt*R) * iStS(i);
        prevCoefXYKALPHA_3 = prevCoefXYKALPHA_3 - im * (iF.*iGv(i,:).*iP);
    end
    coefXYKALPHA_3 = coefXYK_3 * iALPHAcoef;
    
    %System matrix and part of the RHS vector due to previous terms
    A     = (coefXYKALPHA_1 + coefXYKALPHA_2)*matrices(idim).M + coefXYKALPHA_3*matrices(idim).varM;
    prevA = matrices(idim).M * (pgdRB{idim} * (prevCoefXYKALPHA_1 + prevCoefXYKALPHA_2).') + ...
            matrices(idim).varM * (pgdRB{idim} * prevCoefXYKALPHA_3.');
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsK = rhs{dimK}{i};
        rhsXY = rhs{dimXY}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            iF = iF + (Rt*rhsXY{j}) * (St*rhsK{j}).';
        end
        prevCoefXYKALPHA_4 = prevCoefXYKALPHA_4 + icoef*iF;    
    end
    prevA = prevA + iALPHArhscoef*prevCoefXYKALPHA_4*rhs{idim};
    
%------------    
    
elseif mystrcmp({'XY' 'THETA'},restnames)
    
    %Dimensions
    restdimsaux = findPGDdimension({'XY' 'THETA'},parameters.PGDdimensions);
    dimXY       = restdimsaux(1);
    dimTHETA    = restdimsaux(2);
    
    %Some definitions and initializations
    R = mode{dimXY};
    Rt = R';
    T = mode{dimTHETA};
    Tt = T';
    iTt = Tt * matrices(dimTHETA).M;
    iH = iTt * pgdRB{dimTHETA};
    alphaT = Tt * rhs{dimTHETA};
    nOfprevTerms = size(pgdRB{1},2);
    coefXY_1 = 0; coefXYALPHA_2 = 0; coefXY_3 = 0;
    prevCoefXYTALPHA_1 = 0; prevCoefXYTALPHA_2 = 0; prevCoefXYTALPHA_3 = 0;
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
        P      = mode{idimALPHA};
        Pt     = P';
        iPt    = Pt * matrices(idimALPHA).M;
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
            iRt = Rt * matrices(dimXY).(ifields{IC}){i};
            iF = iRt * pgdRB{dimXY};
            coefXY_1 = coefXY_1 + icoefs(IC) * (iRt*R);
            prevCoefXYTALPHA_1 = prevCoefXYTALPHA_1 - icoefs(IC) * (iF.*iPiH);
        end
    end
    coefXYALPHA_1 = coefXY_1 * iALPHAcoef;
    
    %Coefficients for the weak form terms with ALPHA boundaries
    ivarALPHAcoef = ivarPtcoef;
    ivarP = VivarPtcoef;
    for iALPHA = 1:nALPHA
        for jALPHA = setdiff(1:nALPHA,iALPHA)
            ivarALPHAcoef(iALPHA) = ivarALPHAcoef(iALPHA) * iPtcoef(jALPHA);
            ivarP(iALPHA,:) = ivarP(iALPHA,:) .* ViPtcoef(jALPHA,:);
        end
    end
    
    %Add ALPHA boundaries which are not idim (using the last coefficient)
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            iRt = Rt * matrices(dimXY).C_alpha{idimALPHA}{i};
            iF  = iRt * pgdRB{dimXY};
            coefXYALPHA_2 = coefXYALPHA_2 + im * ivarALPHAcoef(iALPHA) * (iRt*R);
            prevCoefXYTALPHA_2 = prevCoefXYTALPHA_2 - im * (iF.*ivarP(iALPHA,:).*iH);
        end
    end
    
    %Coefficient for ALPHA boundary corresponding to idim (using the last coefficient)
    for i = 1:coefs(coefsINDEX(end)).nOfTerms
        iRt = Rt * matrices(dimXY).C_alpha{idim}{i};
        iF  = iRt * pgdRB{dimXY};
        coefXY_3 = coefXY_3 + im * (iRt*R);
        prevCoefXYTALPHA_3 = prevCoefXYTALPHA_3 - im * (iF.*iPiH);
    end
    coefXYALPHA_3 = coefXY_3 * iALPHAcoef;
    
    %System matrix and part of the RHS vector due to previous terms
    coefTHETA = iTt * T;
    A     = coefTHETA*(coefXYALPHA_1 + coefXYALPHA_2)*matrices(idim).M + ...
            coefTHETA*coefXYALPHA_3*matrices(idim).varM;
    prevA = matrices(idim).M * (pgdRB{idim} * (prevCoefXYTALPHA_1 + prevCoefXYTALPHA_2).') + ...
            matrices(idim).varM * (pgdRB{idim} * prevCoefXYTALPHA_3.');
    
    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    iF = 0;
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        icoef = rhscoefs(i);
        jF = 0;
        for j = 1:numel(rhsXY)
            jF = jF + rhsXY{j};
        end
        iF = iF + icoef*jF;  
    end
    prevCoefXYTALPHA_4 = iALPHArhscoef * (Rt*iF) * alphaT.';
    prevA = prevA + prevCoefXYTALPHA_4 * rhs{idim};
    
%------------
    
elseif mystrcmp('XY',restnames)
    
    %Dimensions
    dimXY = findPGDdimension('XY',parameters.PGDdimensions);
    
    %Some definitions and initializations
    R = mode{dimXY};
    Rt = R';
    nOfprevTerms = size(pgdRB{1},2);
    coefXYALPHA_1 = 0; coefXYALPHA_2 = 0; coefXYALPHA_3 = 0;
    prevCoefXYALPHA_1 = 0; prevCoefXYALPHA_2 = 0; prevCoefXYALPHA_3 = 0;
    prevCoefXYALPHA_4 = 0;
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
        P      = mode{idimALPHA};
        Pt     = P';
        iPt    = Pt * matrices(idimALPHA).M;
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
            iF = iRt * pgdRB{dimXY};
            coefXYALPHA_1 = coefXYALPHA_1 + icoefs(IC) * iALPHAcoef * (iRt*R);
            prevCoefXYALPHA_1 = prevCoefXYALPHA_1 - icoefs(IC) * (iF.*iP);
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
    
    %Add ALPHA boundaries which are not idim (using the last coefficient)
    for iALPHA = 1:nALPHA
        idimALPHA = dimALPHA(iALPHA);
        for i = 1:coefs(coefsINDEX(end)).nOfTerms
            iRt = Rt * matrices(dimXY).C_alpha{idimALPHA}{i};
            iF  = iRt * pgdRB{dimXY};
            coefXYALPHA_2 = coefXYALPHA_2 + im * ivarALPHAcoef(iALPHA) * (iRt*R);
            prevCoefXYALPHA_2 = prevCoefXYALPHA_2 - im * (iF.*ivarP(iALPHA,:));
        end
    end
    
    %Coefficient for ALPHA boundary corresponding to idim (using the last coefficient)
    for i = 1:coefs(coefsINDEX(end)).nOfTerms
        iRt = Rt * matrices(dimXY).C_alpha{idim}{i};
        iF  = iRt * pgdRB{dimXY};
        coefXYALPHA_3 = coefXYALPHA_3 + im * iALPHAcoef * (iRt*R);
        prevCoefXYALPHA_3 = prevCoefXYALPHA_3 - im * (iF.*iP);
    end
    
    %System matrix and part of the RHS vector due to previous terms
    A     = (coefXYALPHA_1 + coefXYALPHA_2)*matrices(idim).M + coefXYALPHA_3*matrices(idim).varM;
    prevA = matrices(idim).M * (pgdRB{idim} * (prevCoefXYALPHA_1 + prevCoefXYALPHA_2).') + ...
            matrices(idim).varM * (pgdRB{idim} * prevCoefXYALPHA_3.');

    %Rest of the RHS vector
    rhscoefs = [-1,-1,1,1,1,im];
    for i = 1:numel(rhscoefs)
        rhsXY = rhs{dimXY}{i};
        icoef = rhscoefs(i);
        iF = 0;
        for j = 1:numel(rhsXY)
            iF = iF + rhsXY{j};
        end
        prevCoefXYALPHA_4 = prevCoefXYALPHA_4 + icoef*sum(Rt*iF);    
    end
    prevA = prevA + iALPHArhscoef*prevCoefXYALPHA_4*rhs{idim};
end

%Compute mode P(alpha)
[lu_L,lu_U,lu_P,lu_Q] = lu(A); %Call UMFPACK
P = lu_Q*(lu_U\(lu_L\(lu_P*(prevA))));

%Normalize if needed
if idim < parameters.nOfPGDdimensions
    P2 = P' * matrices(idim).M * P;
    P = P/sqrt(P2);
end


