function coefs = computeCoefs_XY(parameters)

%Some variables
omega = parameters.fixedParametricDims(1);
dimXY = findPGDdimension('XY',parameters.PGDdimensions);

%Data nedeed for coefficients
k         = computeWaveNumber(omega,parameters.meshes(dimXY).bottom);
ccg       = celerities(omega,k,parameters.meshes(dimXY).bottom);
sx        = 1 + 1i * parameters.meshes(dimXY).PML.sigma(:,1) / omega;
sy        = 1 + 1i * parameters.meshes(dimXY).PML.sigma(:,2) / omega;

%Evaluation of coefficients
coefsHandle = {
               @(x,y,sx,sy)x                 %ccg
               @(x,y,sx,sy)x.*(sy./sx)       %ccg * (sigma_y/sigma_x)
               @(x,y,sx,sy)x.*(sx./sy)       %ccg * (sigma_x/sigma_y)
               @(x,y,sx,sy)x.*sx.*sy.*(y.^2) %ccg * sigma_x * sigma_y * k^2;
               @(x,y,sx,sy)x.*y              %ccg * k
              };
nOfCoefs = length(coefsHandle);
coefs = struct();
for i = 1:nOfCoefs 
    coefs(i).RB{dimXY} = coefsHandle{i}(ccg,k,sx,sy); 
    coefs(i).nOfTerms = 1;
end











