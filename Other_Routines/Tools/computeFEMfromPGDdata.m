function ufem = computeFEMfromPGDdata(snapvalue,meshes,parameters)

%PGD dimensions
dimXY = findPGDdimension('XY',parameters.PGDdimensions);

%Snapshot value 
omega   = snapvalue(1);
theta   = snapvalue(2);

%Data nedeed for coefficients
k         = computeWaveNumber(omega,parameters.meshes(dimXY).bottom);
ccg       = celerities(omega,k,parameters.meshes(dimXY).bottom);
sx        = 1 + 1i * parameters.meshes(dimXY).PML.sigma(:,1) / omega;
sy        = 1 + 1i * parameters.meshes(dimXY).PML.sigma(:,2) / omega;

%Evaluation of coefficients at snapshot value
coefsHandle = {
               @(x,y,sx,sy)x                 %ccg
               @(x,y,sx,sy)x.*(sy./sx)       %ccg * (sigma_y/sigma_x)
               @(x,y,sx,sy)x.*(sx./sy)       %ccg * (sigma_x/sigma_y)
               @(x,y,sx,sy)x.*sx.*sy.*(y.^2) %ccg * sigma_x * sigma_y * k^2;
               @(x,y,sx,sy)x.*y              %ccg * k
              };
nOfCoefs = length(coefsHandle);
coefs = cell(nOfCoefs,1);
for i = 1:nOfCoefs, coefs{i} = coefsHandle{i}(ccg,k,sx,sy); end

%Evaluation of incident wave -planar wave- at snapshot value
x    = meshes(dimXY).X(:,1);
y    = meshes(dimXY).X(:,2);
node = meshes(dimXY).T.ext(1);
k0   = k(node);
ip   = exp(1i*k0.*(x.*cos(theta) + y.*sin(theta))); 

%Volume matrices K (Kint, Kxpml, Kypml) and varM -mass matrix with coefficients- (varMint, varMpml)
disp('      Volume matrices...')
[Kint,varMint,Kxpml,Kypml,varMpml] = PGDberkhoffVolumeMatrices(...
    meshes(dimXY).X,...
    meshes(dimXY).T.all,...
    meshes(dimXY).referenceElement,...
    coefs{1},... %for Kint
    coefs{2},... %for Kxpml
    coefs{3},... %for Kypml
    coefs{4},... %for varMint, varMpml
    parameters.meshes(dimXY).PML.elements);
varM{1} = varMint{1} + varMpml{1};

%Boundary mass matrix alpha*a(x,y)*Ni*Nj in gamma_R and gamma_PML with alpha = const for each element
disp('      Boundary matrices...')
if ~isempty(meshes(dimXY).Tb.gammaR)
    Ccell = PGDmassMatrix1D(...
        meshes(dimXY).X,...
        meshes(dimXY).Tb.gammaR,...
        meshes(dimXY).referenceElement,...
        {coefs(5)},...
        dimXY,... %XY dimension index of coefficients, i.e. RB{1}
        parameters.meshes(dimXY).alpha);
    C = Ccell{1}{1};
else
    sizeC = size(meshes(dimXY).X,1);
    C = spalloc(sizeC,sizeC,0);
end
namesALPHA = fieldnames(meshes(dimXY).Tb.gammaALPHA);
nALPHA = numel(namesALPHA);
for j = 1:nALPHA
    jparam = snapvalue(2+j); %alpha value of the snapshot
    Ccell_alpha = PGDmassMatrix1D(...
        meshes(dimXY).X,...
        meshes(dimXY).Tb.gammaALPHA.(namesALPHA{j}),...
        meshes(dimXY).referenceElement,...
        {coefs(5)},...
        dimXY,... %XY dimension index of coefficients, i.e. RB{1}
        jparam*ones(size(meshes(dimXY).Tb.gammaPML,1),1));
    C = C + Ccell_alpha{1}{1};
end
Dcell = PGDmassMatrix1D(...
    meshes(dimXY).X,...
    meshes(dimXY).Tb.gammaPML,...
    meshes(dimXY).referenceElement,...
    {coefs(5)},...
    dimXY,...
    ones(size(meshes(dimXY).Tb.gammaPML,1),1)); %alpha = 1
D = Dcell{1};

%RHS vectors
disp('      Right hand side vectors...')

%Vectors in \Omega_PML
rhs1 = Kxpml{1}   * ip;
rhs2 = Kypml{1}   * ip;
rhs3 = varMpml{1} * ip;

%Vectors in \partial\Gamma_rest
[~,irow] = intersect(meshes(dimXY).Tb.boundaryPML.conec,meshes(dimXY).Tb.gammaPML,'rows');
elementsPML_rest = setdiff(1:size(meshes(dimXY).Tb.boundaryPML.conec,1),irow);
[rhs4,rhs5] = PGDberkhoffBoundaryVectors(...
    meshes(dimXY).X,...
    meshes(dimXY).T.all,...
    meshes(dimXY).Tb.boundaryPML.conec(elementsPML_rest,:),...
    meshes(dimXY).referenceElement,...
    meshes(dimXY).elementFaceInfo.boundaryPML(elementsPML_rest,1),...
    meshes(dimXY).Tb.boundaryPML.der2DShapeFunOn1D(:,:,:,elementsPML_rest),...
    ip,...
    coefs{2},...
    coefs{3});

%Vectors in \Gamma_PML
rhs6 = D{1} * ip;

%FEM solution at snapshot value
disp('      Solving the linear system...')
sysMat    = -(Kxpml{1} + Kypml{1} + Kint{1}) + varM{1} + 1i*(C + D{1});
sysVec    = -(rhs1 + rhs2) + rhs3 + rhs4{1} + rhs5{1} + 1i*rhs6;
[L,U,P,Q] = lu(sysMat); %Call UMFPACK
ufem      = Q*(U\(L\(P*sysVec)));





