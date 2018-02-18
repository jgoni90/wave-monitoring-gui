
%% General

clear all
home
parameters.PGDdimensions    = {'K' ; 'DEPTH'};
parameters.nOfPGDdimensions = 2;

%% Frequency / depth data

aT              = 3;
bT              = 16;
ad              = 1;
bd              = 65;
ininOfNodes_w   = 1000;
nDeg_w          = 1;
ininOfNodes_d   = 1000;
nDeg_d          = 1;

%% Data for PGD projection

OnlyPlotTensor                      = 0;

runSVD                              = false;
parameters.maxterms                 = 100;
algorithm.projection.residualType   = 'TENSOR_APP';
parameters.residualEachTerm         = 1;
algorithm.projection.maxiter        = 20;
parameters.toliter                  = 1e-3;
algorithm.projection.relerror       = 1e-15;
algorithm.projection.projCoeff      = true;
parameters.recoverPGD.value         = 0;
parameters.recoverPGD.file          = '';
parameters.outputFileName           = 'PGD_freq_depth_OUTPUT.mat';

%% Meshes dimension FREQUENCY x DEPTH

[meshes(1).X,meshes(1).T] = createPGDMesh('K',...
    struct('Tini',aT,'Tend',bT,'nDeg',nDeg_w,'ininOfNodes',ininOfNodes_w));
nOfNodes_w = size(meshes(1).X,1);
meshes(1).referenceElement = createReferenceElement(1,(nDeg_w+1)*(nDeg_w+2)/2);

[meshes(2).X,meshes(2).T] = createPGDMesh('THETA',...
    struct('THETAini',ad,'THETAend',bd,'nDeg',nDeg_d,'ininOfNodes',ininOfNodes_d));
nOfNodes_d = size(meshes(2).X,1);
meshes(2).referenceElement = createReferenceElement(1,(nDeg_d+1)*(nDeg_d+2)/2);

%% Mass matrices

for i = 1:2 
    auxones = {ones(size(meshes(i).X,1),1)};
    Mw = PGDmassMatrix1D(...
        meshes(i).X,...
        meshes(i).T,...
        meshes(i).referenceElement,...
        {auxones},...   
        1,...
        ones(size(meshes(i).T,1),1));
    matrices(i).M = Mw{1}{1};
end

%% Data tensors and PGD projection

nOfTensors = 3;
coefsHandle = {@(x,y)x, @(x,y)x.*(y.^2), @(x,y)x.*y};

for i = 1:nOfTensors
    
    %Data tensor
    t = zeros(nOfNodes_w,nOfNodes_d);
    for j = 1:nOfNodes_w
        k = computeWaveNumber(meshes(1).X(j),meshes(2).X);
        ccg = celerities(meshes(1).X(j),k,meshes(2).X);
        t(j,:) = coefsHandle{i}(ccg,k);
    end
    
    %Plots
    if OnlyPlotTensor
        
        %Surf plot
        figure
        [XX,YY] = meshgrid(meshes(2).X,meshes(1).X);
        hs = surf(XX,YY,t);
        set(hs,'edgealpha',0)
        set(gca,'fontsize',30)
        xlim([min(meshes(2).X),max(meshes(2).X)])
        ylim([min(meshes(1).X),max(meshes(1).X)])
        ylabel('frequency'), xlabel('bottom depth')
        box on
        view(2)
        hold on
        
        %Contour lines plot
        [Cplot,hc] = contour3(XX,YY,t,'k-');
        set(hc,'linewidth',3)
        
        %Select contour labels
        ht = clabel(Cplot,hc,'manual');
        set(ht,'fontsize',20)
        
        continue
    end
    
    %Shape functions tensor for PGD
    A = tensor_tijNiNj(meshes(1).X,meshes(1).T,meshes(2).X,meshes(2).T,...
        meshes(1).referenceElement,meshes(2).referenceElement,t);
    
    %PGD projection
    printPGD()
    pgd(i) = PGDprojection_NonSeparableTensor(...
                A,...
                matrices,...
                meshes,...
                t,...
                parameters,...
                algorithm);
            
    %Perform SVD
    if runSVD
        disp('Running SVD...')
        [U,s,V] = svd(t,0); %economic size SVD
        sdiag = diag(s);
        SVDalpha{i} = sdiag;
    else
        SVDalpha{i} = [];
    end
end

%% Save

if ~OnlyPlotTensor, save(parameters.outputFileName,'pgd','parameters','meshes','matrices','SVDalpha'), end










