
%% General

%First manually load PGDmeshes and parameters
%This routine evaluates the PGDprojection in the entire mesh (not optimal)

%clear all
home
dimXY       = findPGDdimension('XY',parameters.PGDdimensions);
dimK        = findPGDdimension('K',parameters.PGDdimensions);
nOfNodes_XY = size(PGDmeshes(dimXY).X,1);
nOfNodes_K  = size(PGDmeshes(dimK).X,1);
ones1       = ones(size(PGDmeshes(1).X,1),1);

%% Data for PGD projection

parametersaux.nOfPGDdimensions          = 2;
parametersaux.maxterms                  = 100;
algorithm.projection.residualType       = 'TENSOR_APP';
parametersaux.residualEachTerm          = 1;
algorithm.projection.maxiter            = 20;
parametersaux.toliter                   = 1e-3;
algorithm.projection.relerror           = 1e-12;
algorithm.projection.projCoeff          = true;
parametersaux.recoverPGD.value          = 0;
parametersaux.recoverPGD.file           = '';
parametersaux.outputFileName            = 'PGD_freq_depth_OUTPUT.mat';

%% Mass matrices

auxones = ones(size(PGDmeshes(dimXY).X,1),1);
[~,Mint,~,~,Mpml] = PGDberkhoffVolumeMatrices(...
    PGDmeshes(dimXY).X,...
    PGDmeshes(dimXY).T.all,...
    PGDmeshes(dimXY).referenceElement,...
    auxones,...
    auxones,...
    auxones,...
    auxones,...
    []); %no PML elements
matrices(dimXY).M = Mint{1} + Mpml{1};

auxones = {ones(size(PGDmeshes(dimK).X,1),1)};
Mw = PGDmassMatrix1D(...
    PGDmeshes(dimK).X,...
    PGDmeshes(dimK).T,...
    PGDmeshes(dimK).referenceElement,...
    {auxones},...   
    1,...
    ones(size(PGDmeshes(dimK).T,1),1));
matrices(dimK).M = Mw{1}{1};

%% Data tensors and PGD projection

%Tensors: analytical definition
coefsHandle = {
%                @(x,y,sx,sy,g)x                 %ccg
%                @(x,y,sx,sy,g)x.*(sy./sx)       %ccg * (sigma_y/sigma_x)
%                @(x,y,sx,sy,g)x.*(sx./sy)       %ccg * (sigma_x/sigma_y)
%                @(x,y,sx,sy,g)x.*sx.*sy.*(y.^2) %ccg * sigma_x * sigma_y * k^2;
%                @(x,y,sx,sy,g)x.*y              %ccg * k
%                @(x,y,sx,sy,g)x.*y.*g           %ccg * k * gamma
               @(x,y,k,t)exp(1i*k.*(x.*cos(t) + y.*sin(t)));   %incident wave
              };
nOfTensors  = length(coefsHandle);

%Tensors: spatial domain
coefsXYfield = {
                %parameters
%                 'T'  'int'                 
%                 'T'  'ext'
%                 'T'  'ext'
%                 'T'  'all'
%                 'Tb' 'gammaR'
%                 'Tb' 'gammaPML'
                
                %incident wave
                'T'  'ext'
               };
          
%Evaluated tensors
disp('Evaluating tensors...')
t = zeros(nOfNodes_K,nOfNodes_XY,nOfTensors);   
for i = 1:nOfTensors
%     xynodes = unique(PGDmeshes(dimXY).(coefsXYfield{i,1}).(coefsXYfield{i,2}));
    xynodes = 1:nOfNodes_XY;
    for j = 1:nOfNodes_K
        
        %Parameters
%         k = computeWaveNumber(PGDmeshes(dimK).X(j),parameters.meshes(dimXY).bottom(xynodes));
%         ccg = celerities(PGDmeshes(dimK).X(j),k,parameters.meshes(dimXY).bottom(xynodes));
%         sx = 1 + 1i * parameters.meshes(dimXY).PML.sigma(xynodes,1) / PGDmeshes(dimK).X(j);
%         sy = 1 + 1i * parameters.meshes(dimXY).PML.sigma(xynodes,2) / PGDmeshes(dimK).X(j);
%         g = 1 + 1i * parameters.meshes(dimXY).sigmaOfGammaCoef(xynodes) / PGDmeshes(dimK).X(j);
%         t(j,xynodes,i) = coefsHandle{i}(ccg,k,sx,sy,g);
        
        %Incident wave
        k = computeWaveNumber(PGDmeshes(dimK).X(j),parameters.meshes(dimXY).bottom(xynodes));
        x = PGDmeshes(dimXY).X(xynodes,1);
        y = PGDmeshes(dimXY).X(xynodes,2);
        theta = parameters.fixedParametricDims(2);
        t(j,xynodes,i) = coefsHandle{i}(x,y,k,theta);
    end         
end

%Shape functions tensor for PGD
disp('Building tensor integrals...')
A = tensor2_tijkNiNjk(PGDmeshes(dimK).X,PGDmeshes(dimK).T,...
    PGDmeshes(dimXY).X,PGDmeshes(dimXY).T.all,...
    PGDmeshes(dimK).referenceElement,PGDmeshes(dimXY).referenceElement,t);

%PGD projection
for i = 1:nOfTensors
    printPGD()
    pgd(i) = PGDprojection_NonSeparableTensor(...
                A(:,:,i).',...
                matrices,...
                PGDmeshes,...
                t(:,:,i).',...
                parametersaux,...
                algorithm);
end

%Add number of terms
for i = 1:nOfTensors
    pgd(i).nOfTerms = size(pgd(i).RB{1},2);
end










