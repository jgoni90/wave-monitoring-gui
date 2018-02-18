function [matrices,rhs] = computePGDMatrices(meshes,coefs,parameters,algorithm,IW)

%% COMPUTE MATRICES FOR EACH DIMENSION

dimX        = findPGDdimension('X',parameters.PGDdimensions);
dimXY       = findPGDdimension('XY',parameters.PGDdimensions);
dimK        = findPGDdimension('K',parameters.PGDdimensions);
dimTHETA    = findPGDdimension('THETA',parameters.PGDdimensions);
dimKTHETA   = findPGDdimension('KTHETA',parameters.PGDdimensions);
dimALPHA    = findPGDdimension('ALPHA*',parameters.PGDdimensions);

matrices = struct();
rhs = cell(1);

%------

if dimX
    
end

%------

if dimXY
    
    disp('  Dimension XY')
    
    %Volume matrices K (Kint, Kxpml, Kypml) and varM -mass matrix with coefficients- (varMint, varMpml)
    disp('      Volume matrices...')
    [Kint,varMint,Kxpml,Kypml,varMpml] = PGDberkhoffVolumeMatrices(...
        meshes(dimXY).X,...
        meshes(dimXY).T.all,...
        meshes(dimXY).referenceElement,...
        coefs(1).RB{dimXY},... %for Kint
        coefs(2).RB{dimXY},... %for Kxpml
        coefs(3).RB{dimXY},... %for Kypml
        coefs(4).RB{dimXY},... %for varMint, varMpml
        parameters.meshes(dimXY).PML.elements);
    varM = cell(coefs(4).nOfTerms);
    for j = 1:coefs(4).nOfTerms, varM{j} = varMint{j} + varMpml{j}; end
    
    %Mass matrix
    auxones = ones(size(meshes(dimXY).X,1),1);
    [~,Mint,~,~,Mpml] = PGDberkhoffVolumeMatrices(...
        meshes(dimXY).X,...
        meshes(dimXY).T.all,...
        meshes(dimXY).referenceElement,...
        auxones,...
        auxones,...
        auxones,...
        auxones,...
        parameters.meshes(dimXY).PML.elements);
    M = Mint{1} + Mpml{1};
    
    %Boundary mass matrix alpha*a(x,y)*Ni*Nj in gamma_ALPHA, gamma_R and gamma_PML with alpha = const for each element
    disp('      Boundary matrices...')
    
    %For variable reflecting coefficient ALPHA
    namesALPHA = fieldnames(meshes(1).Tb.gammaALPHA);
    nALPHA = numel(namesALPHA);
    C_alpha = cell(parameters.nOfPGDdimensions,1);
    for j = 1:nALPHA
        Ccell_alpha = PGDmassMatrix1D(...
            meshes(dimXY).X,...
            meshes(dimXY).Tb.gammaALPHA.(namesALPHA{j}),...
            meshes(dimXY).referenceElement,...
            {coefs(5).RB},...
            dimXY,... %XY dimension index of coefficients, i.e. RB{1}
            ones(size(meshes(dimXY).Tb.gammaPML,1),1)); %alpha = 1
        C_alpha{dimALPHA(j)} = Ccell_alpha{1};
    end
    
    %For fixed reflecting coefficient ALPHA (boundary gammaR) 
    if ~isempty(meshes(dimXY).Tb.gammaR)
        Ccell = PGDmassMatrix1D(...
            meshes(dimXY).X,...
            meshes(dimXY).Tb.gammaR,...
            meshes(dimXY).referenceElement,...
            {coefs(5).RB},...
            dimXY,... %XY dimension index of coefficients, i.e. RB{1}
            parameters.meshes(dimXY).alpha);
    else
        Ccell = {[]};
    end
    C = Ccell{1};
    
    %For PML boundary gammaPML
    Dcell = PGDmassMatrix1D(...
        meshes(dimXY).X,...
        meshes(dimXY).Tb.gammaPML,...
        meshes(dimXY).referenceElement,...
        {coefs(5).RB},...
        dimXY,...
        ones(size(meshes(dimXY).Tb.gammaPML,1),1)); %alpha = 1
    D = Dcell{1};
    
    %RHS vectors
    disp('      Right hand side vectors...')
    
    %Sparse extension of the XY reduced basis on the PML to the entire XY mesh 
    f1 = zeros(parameters.nOfPGDdimensions,1);
    f1(dimXY) = true;
    f2 = cell(parameters.nOfPGDdimensions,1);
    f2{dimXY}{1} = IW.meshes(dimXY).nodes;
    f2{dimXY}{2} = size(meshes(dimXY).X,1);
    ip = arrangePGD(IW,[],f1,f2,'extends');
    
    %Initialize
    nOfRHSterms = 6;
    rhs{dimXY} = cell(nOfRHSterms,1);
    auxv = [2,3,4,2,3,5];
    for i = 1:nOfRHSterms, rhs{dimXY}{i} = cell(coefs(auxv(i)).nOfTerms,1); end
        
    %Vectors in \Omega_PML
    for j = 1:coefs(2).nOfTerms, rhs{dimXY}{1}{j} = Kxpml{j}    * ip.RB{dimXY}; end
    for j = 1:coefs(3).nOfTerms, rhs{dimXY}{2}{j} = Kypml{j}    * ip.RB{dimXY}; end
    for j = 1:coefs(4).nOfTerms, rhs{dimXY}{3}{j} = varMpml{j}  * ip.RB{dimXY}; end

    %Vectors in \partial\Gamma_rest
    [~,irow] = intersect(meshes(dimXY).Tb.boundaryPML.conec,meshes(dimXY).Tb.gammaPML,'rows');
    elementsPML_rest = setdiff(1:size(meshes(dimXY).Tb.boundaryPML.conec,1),irow);
    [rhs{dimXY}{4},rhs{dimXY}{5}] = PGDberkhoffBoundaryVectors(...
        meshes(dimXY).X,...
        meshes(dimXY).T.all,...
        meshes(dimXY).Tb.boundaryPML.conec(elementsPML_rest,:),...
        meshes(dimXY).referenceElement,...
        meshes(dimXY).elementFaceInfo.boundaryPML(elementsPML_rest,1),...
        meshes(dimXY).Tb.boundaryPML.der2DShapeFunOn1D(:,:,:,elementsPML_rest),...
        ip.RB{dimXY},...
        coefs(2).RB{dimXY},...
        coefs(3).RB{dimXY});
    
    %Vectors in \Gamma_PML
    for j = 1:coefs(5).nOfTerms, rhs{dimXY}{6}{j} = D{j} * ip.RB{dimXY}; end
    
    clear ip
    
    %Store
    matrices(dimXY).Kint    = Kint;
    matrices(dimXY).Kxpml   = Kxpml;
    matrices(dimXY).Kypml   = Kypml;
    matrices(dimXY).varM    = varM;
    matrices(dimXY).M       = M;
    matrices(dimXY).C_alpha = C_alpha;
    matrices(dimXY).C       = C;
    matrices(dimXY).D       = D;
    
    clear Kint Kxpml Kypml varM M C_alpha C D
end

%------

if dimK
    
    disp('  Dimension OMEGA')
    
    %Volume matrices varM -mass matrix with coefficients-
    disp('      Volume matrices...')
    varMw = PGDmassMatrix1D(...
        meshes(dimK).X,...
        meshes(dimK).T,...
        meshes(dimK).referenceElement,...
        {coefs.RB},...  %Mass matrix for all coefficients
        dimK,...        %Omega dimension index of coefficients, i.e. RB{2}
        ones(size(meshes(dimK).T,1),1));
    
    %Mass matrix
    auxones = {ones(size(meshes(dimK).X,1),1)};
    Mw = PGDmassMatrix1D(...
        meshes(dimK).X,...
        meshes(dimK).T,...
        meshes(dimK).referenceElement,...
        {auxones},...
        1,...
        ones(size(meshes(dimK).T,1),1));
    
    %RHS vectors
    nOfRHSterms = 6;
    rhs{dimK} = cell(nOfRHSterms,1);
    auxv = [2,3,4,2,3,5];
    if dimTHETA %in \Omega_{omega x theta}
        disp('      Right hand side vectors in dimension OMEGA x THETA...')
        
        nOfNodesK      = size(meshes(dimK).X,1);
        nOfNodesTHETA  = size(meshes(dimTHETA).X,1);
        nOfNodesKTHETA = size(IW.meshes(dimK).X,1);
        auxonesTHETA   = ones(1,nOfNodesTHETA,1);
        auxonesKTHETA  = ones(nOfNodesKTHETA,1);
        for i = 1:nOfRHSterms
            
            %Matrices in QUAD mesh KTHETA for each coefficient
            icoef = auxv(i);
            nOfcoefTerms = coefs(icoef).nOfTerms;
            icoefKTHETA = bsxfun(@times,reshape(coefs(icoef).RB{dimK},nOfNodesK,1,nOfcoefTerms),...
                auxonesTHETA);
            [~,varMwt] = PGDberkhoffVolumeMatrices(...
                IW.meshes(dimK).X,...
                IW.meshes(dimK).T,...
                IW.meshes(dimK).referenceElement,...
                auxonesKTHETA,...
                auxonesKTHETA,...
                auxonesKTHETA,...
                reshape(icoefKTHETA,nOfNodesK*nOfNodesTHETA,nOfcoefTerms),... %for varMwt
                []); %no PML elements
            
            %Vectors in dimension KTHETA for each coefficient and all IW terms
            rhs{dimK}{i} = cell(nOfcoefTerms,1);
            for j = 1:nOfcoefTerms, rhs{dimK}{i}{j} = varMwt{j} * IW.RB{dimK}; end
        end
        
        clear varMwt
        
    else %in \Omega_omega
        disp('      Right hand side vectors...')
        
        for i = 1:nOfRHSterms
            icoef = auxv(i);
            nOfcoefTerms = coefs(icoef).nOfTerms;
            rhs{dimK}{i} = cell(nOfcoefTerms,1);
            for j = 1:nOfcoefTerms, rhs{dimK}{i}{j} = varMw{icoef}{j} * IW.RB{dimK}; end
        end
    end
    
    %Store
    matrices(dimK).varM = varMw;
    matrices(dimK).M    = Mw{1}{1};
    
    clear varMw Mw
end

%------

if dimKTHETA
    
    disp('  Dimension OMEGA x THETA')
    
    %Volume matrices varM -mass matrix with coefficients-
    disp('      Volume matrices...')
    nOfNodesK = parameters.meshes(dimKTHETA).nK + 1;
    nOfNodesTHETA = parameters.meshes(dimKTHETA).nTHETA + 1;
    nOfNodesKTHETA = size(meshes(dimKTHETA).X,1);
    auxonesTHETA  = ones(1,nOfNodesTHETA,1);
    auxonesKTHETA = ones(nOfNodesKTHETA,1);
    varMw = cell(1,numel(coefs));
    for i = 1:numel(coefs)
        nOfcoefTerms = coefs(i).nOfTerms;
        icoefKTHETA = bsxfun(@times,reshape(coefs(i).RB{dimKTHETA},nOfNodesK,1,nOfcoefTerms),...
                auxonesTHETA);
        [~,varMw{i}] = PGDberkhoffVolumeMatrices(...
            IW.meshes(dimKTHETA).X,...
            IW.meshes(dimKTHETA).T,...
            IW.meshes(dimKTHETA).referenceElement,...
            auxonesKTHETA,...
            auxonesKTHETA,...
            auxonesKTHETA,...
            reshape(icoefKTHETA,nOfNodesKTHETA,nOfcoefTerms),... %for varM
            []); %no PML elements
    end
    
    %RHS vectors in \Omega_KTHETA
    disp('      Right hand side vectors...')
    nOfRHSterms = 6;
    rhs{dimKTHETA} = cell(nOfRHSterms,1);
    auxv = [2,3,4,2,3,5];
    for i = 1:nOfRHSterms
        icoef = auxv(i);
        nOfcoefTerms = coefs(icoef).nOfTerms;
        rhs{dimKTHETA}{i} = cell(nOfcoefTerms,1);
        for j = 1:nOfcoefTerms, rhs{dimKTHETA}{i}{j} = varMw{icoef}{j} * IW.RB{dimKTHETA}; end
    end
    
    %Store
    matrices(dimKTHETA).varM = varMw;
    matrices(dimKTHETA).M    = IW.matrices(dimKTHETA).M;
    
    clear varMw
end

%------

if dimTHETA
    
    disp('  Dimension THETA')
    
    %Mass matrix
    disp('      Volume matrices...')
    auxones = {ones(size(meshes(dimTHETA).X,1),1)};
    Mt = PGDmassMatrix1D(...
        meshes(dimTHETA).X,...
        meshes(dimTHETA).T,...
        meshes(dimTHETA).referenceElement,...
        {auxones},...   
        1,...
        ones(size(meshes(dimTHETA).T,1),1));
    
    %RHS vectors in \Omega_theta
    if ~dimK
        disp('      Right hand side vectors...')
        
        nOfcoefTerms = size(IW.RB{2},2);
        rhs{dimTHETA} = zeros(size(meshes(dimTHETA).X,1),nOfcoefTerms);
        varMt = PGDmassMatrix1D(...
            meshes(dimTHETA).X,...
            meshes(dimTHETA).T,...
            meshes(dimTHETA).referenceElement,...
            {IW.RB},...   
            dimTHETA,...
            ones(size(meshes(dimTHETA).T,1),1));
        for j = 1:nOfcoefTerms, rhs{dimTHETA}(:,j) = varMt{1}{j} * auxones{1}; end
        
        clear varMt
    end
    
    %Store
    matrices(dimTHETA).M = Mt{1}{1};
    
    clear Mt
end

%------

if dimALPHA %multiple ALPHA dimensions allowed
    
    for idimALPHA = dimALPHA %For each ALPHA dimension
        disp(['  Dimension ' parameters.PGDdimensions{idimALPHA,1}])

        %Volume matrices varM -mass matrix with coefficients-
        disp('      Volume matrices...')
        alphacoef = {meshes(idimALPHA).X};
        varM_alpha = PGDmassMatrix1D(...
            meshes(idimALPHA).X,...
            meshes(idimALPHA).T,...
            meshes(idimALPHA).referenceElement,...
            {alphacoef},...   
            1,...
            ones(size(meshes(idimALPHA).T,1),1));
        
        %Mass matrix
        auxones = {ones(size(meshes(idimALPHA).X,1),1)};
        M_alpha = PGDmassMatrix1D(...
            meshes(idimALPHA).X,...
            meshes(idimALPHA).T,...
            meshes(idimALPHA).referenceElement,...
            {auxones},...
            1,...
            ones(size(meshes(idimALPHA).T,1),1));

        %%RHS vectors in \Omega_alpha
        disp('      Right hand side vectors...')
        rhs{idimALPHA} = M_alpha{1}{1} * auxones{1};

        %Store
        matrices(idimALPHA).varM = varM_alpha{1}{1};
        matrices(idimALPHA).M = M_alpha{1}{1};
    end
end




    