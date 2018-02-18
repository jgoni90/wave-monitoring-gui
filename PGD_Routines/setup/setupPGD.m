%SCRIPT TO SETUP THE PGD CODE FROM THE PGDdata.m file

%% GENERAL PARAMETERS

parameters.PGDdimensions        = PGDdimensions;
parameters.maxterms             = maxterms;
parameters.maxiter              = maxiter;
parameters.toliter              = toliter;
parameters.tolu                 = tolu;
parameters.residualType         = residualType;
parameters.residualEachTerm     = residualEachTerm;
parameters.saveEachTerm         = saveEachTerm;
parameters.outputFile           = outputFile;
parameters.recoverPGD.value     = recoverPGD;
parameters.recoverPGD.file      = recoverFile;
parameters.fixedParametricDims  = fixedParametricDims;
parameters.iwparam              = iwparam;
parameters.coefparam            = iwparam;

%% MESHES

%Meshes and algorithm
algorithmstr      = 'PGDmapping';
updatestr         = 'PGDupdate';
optionstr{1}      = {'dual' 'dual' 'value'};
optionstr{2}      = {'pg' 'PG'};

parameters.nOfPGDdimensions = size(PGDdimensions,1);
PGDmeshes = struct();
PGDalgorithm.function.mapping = cell(parameters.nOfPGDdimensions,1);
PGDalgorithm.function.updateInfo = cell(parameters.nOfPGDdimensions,1);
PGDalgorithm.update.function = cell(parameters.nOfPGDdimensions,1);
parameters.meshes = struct();

disp('BUILDING PGD MESHES')

for idim = 1:parameters.nOfPGDdimensions
    
    dim = PGDdimensions(idim,:);
    namedim = dim{1};
    loadfile = dim{3};
    
    switch namedim
        
        case 'X'
        
        case 'XY'
            
            disp('  Dimension XY...')
            
            if loadfile
                file = load(fileloadmeshes{idim});
                
                PGDmeshes(idim).X = file.data.mesh.X;
                PGDmeshes(idim).T.all = file.data.mesh.T;
                PGDmeshes(idim).T.ext = file.data.mesh.extT;
                PGDmeshes(idim).T.int = file.data.mesh.intT;
                PGDmeshes(idim).referenceElement = file.data.mesh.referenceElement;
                
                %Attributes of the boundaries & bottom
                parameters.meshes(idim).BC = file.data.BC;
                parameters.meshes(idim).boundaryNames = file.data.mesh.boundaryNames;
                parameters.meshes(idim).nOfBoundaries = length(parameters.meshes(idim).boundaryNames);
                parameters.meshes(idim).bottom = file.data.bottom.value;
                parameters.meshes(idim).alpha = [];
                
                %Boundary & PML attributes
                parameters.meshes(idim).PML.elements = file.data.mesh.inPMLelems;
                PGDmeshes(idim).elementFaceInfo = file.data.mesh.elementFaceInfo;
                PGDmeshes(idim).Tb.boundaryPML.conec = file.data.mesh.extTb;
                PGDmeshes(idim).Tb.boundaryPML.der2DShapeFunOn1D = file.data.mesh.der2DShapeFunOn1D.extTb;
                PGDmeshes(idim).Tb.gammaR = [];
                PGDmeshes(idim).Tb.gammaPML = [];
                PGDmeshes(idim).Tb.gammaALPHA = struct();
                PGDmeshes(idim).Tb.all = [];
                gammaCont = 0;
                fixedParametricDimsAUX = parameters.fixedParametricDims;
                for i = 1:parameters.meshes(idim).nOfBoundaries
                    iname = parameters.meshes(idim).boundaryNames{i};
                    iTb = file.data.mesh.(['Tb_' iname]);
                    xy_nodes = unique(iTb);
                    icond = parameters.meshes(idim).BC.values(i);
                    
                    %Particularize depending on the boundary condition
                    if icond == 4 %PML (assumed only one PML)
                        parameters.meshes(idim).PML.sigma = file.data.PML{1,i};
                        PGDmeshes(idim).Tb.gammaPML = [PGDmeshes(idim).Tb.gammaPML ; iTb];
                    
                    elseif icond == 3 %Radiation (assumed in the PML)
                        PGDmeshes(idim).Tb.gammaPML = [PGDmeshes(idim).Tb.gammaPML ; iTb];
                    
                    elseif icond == 2 && ... %Variable reflecting (assumed outside the PML)
                        mystrcmp(['ALPHA_' iname],parameters.PGDdimensions(:,1))
                        gammaCont = gammaCont + 1;
                        PGDmeshes(idim).Tb.gammaALPHA.(iname) = iTb;
                        
                    elseif icond == 2 %Fixed reflecting (assumed outside the PML)
                        gammaCont = gammaCont + 1;
                        PGDmeshes(idim).Tb.gammaR = [PGDmeshes(idim).Tb.gammaR ; iTb];
                        if parameters.fixedParametricDims(3) == -2
                            iparam = parameters.meshes(idim).BC.parameters{i}{icond};
                            fixedParametricDimsAUX(gammaCont+2) = iparam;
                        else
                            iparam = parameters.fixedParametricDims(gammaCont+2);
                        end
                        parameters.meshes(idim).alpha = [parameters.meshes(idim).alpha ; iparam*ones(size(iTb,1),1)];  
                    end
                    
                    %Entire boundary
                    PGDmeshes(idim).Tb.all = [PGDmeshes(idim).Tb.all ; iTb];
                end
                parameters.fixedParametricDims = fixedParametricDimsAUX;
                
                %Dual connectivity
                if PGDalgorithm.dual.value
                    PGDmeshes(idim).Tdual = file.data.mesh.(PGDalgorithm.dual.spatialAttributeName);
                end
                
                clear file
            else
                error('XY meshes have to be loaded from Berkhoff GUI file')
            end
            
            namedimFile = namedim;
            
        case 'K'
            
            disp('  Dimension OMEGA...')
            
            meshparams = meshparameters{idim};
            
            parameters.meshes(idim).Tini            = meshparams{1};
            parameters.meshes(idim).Tend            = meshparams{2};
            parameters.meshes(idim).ininOfNodes     = meshparams{3};
            parameters.meshes(idim).nDeg            = meshparams{4};
            parameters.meshes(idim).nDegNodes       = meshparams{5};
            
            if loadfile
                PGDmeshes(idim) = load(fileloadmeshes{idim});
            else
                [PGDmeshes(idim).X,PGDmeshes(idim).T] = createPGDMesh('K',parameters.meshes(idim));
                PGDmeshes(idim).referenceElement = createReferenceElement(1,...
                    parameters.meshes(idim).nDegNodes,[]);
            end
            
            namedimFile = namedim;
            
        case 'THETA'
            
            disp('  Dimension THETA...')
            
            meshparams = meshparameters{idim};
            
            parameters.meshes(idim).THETAini        = meshparams{1}*pi/180;
            parameters.meshes(idim).THETAend        = meshparams{2}*pi/180;
            parameters.meshes(idim).ininOfNodes     = meshparams{3};
            parameters.meshes(idim).nDeg            = meshparams{4};
            parameters.meshes(idim).nDegNodes       = meshparams{5};
            
            if loadfile
                PGDmeshes(idim) = load(fileloadmeshes{idim});
            else
                [PGDmeshes(idim).X,PGDmeshes(idim).T] = createPGDMesh('THETA',parameters.meshes(idim));
                PGDmeshes(idim).referenceElement = createReferenceElement(1,...
                    parameters.meshes(idim).nDegNodes,[]);
            end
            
            namedimFile = namedim;
            
        case 'KTHETA'
            
            disp('  Dimension OMEGA x THETA...')
            
            meshparams = meshparameters{idim};

            parameters.meshes(idim).Tini = meshparams{1};
            parameters.meshes(idim).Tend = meshparams{2};
            parameters.meshes(idim).nK = meshparams{3};
            parameters.meshes(idim).THETAini = meshparams{4}*pi/180;
            parameters.meshes(idim).THETAend = meshparams{5}*pi/180;
            parameters.meshes(idim).nTHETA = meshparams{6};
            parameters.meshes(idim).nDeg = meshparams{7};
            parameters.meshes(idim).nDegNodes = meshparams{8};
            parameters.meshes(idim).equalSpaced = meshparams{9};

            if loadfile
                PGDmeshes(idim) = load(fileloadmeshes{idim});
            else
                [PGDmeshes(idim).X,PGDmeshes(idim).T] = createPGDMesh('KTHETA',parameters.meshes(idim));
                PGDmeshes(idim).referenceElement = createReferenceElement(0,...
                    parameters.meshes(idim).nDegNodes,[]);
            end
            
            namedimFile = namedim;
            
        otherwise %ALPHA dimensions

            disp(['  Dimension ' namedim '...'])

            meshparams = meshparameters{idim};

            parameters.meshes(idim).THETAini        = meshparams{1};
            parameters.meshes(idim).THETAend        = meshparams{2};
            parameters.meshes(idim).ininOfNodes     = meshparams{3};
            parameters.meshes(idim).nDeg            = meshparams{4};
            parameters.meshes(idim).nDegNodes       = meshparams{5};
            parameters.meshes(idim).XYboundaryName  = namedim(7:end);

            if loadfile
                PGDmeshes(idim) = load(fileloadmeshes{idim});
            else
                [PGDmeshes(idim).X,PGDmeshes(idim).T] = createPGDMesh('THETA',parameters.meshes(idim));
                PGDmeshes(idim).referenceElement = createReferenceElement(1,...
                    parameters.meshes(idim).nDegNodes,[]);
            end
            
            namedimFile = 'ALPHA';
    end
    
    %PGD algorithm handle
    namefunction = algorithmstr;
    updatenamefunction = updatestr;
    for cont = 1:numel(optionstr)
        val = PGDalgorithm;
        for cont2 = 2:numel(optionstr{cont})
            val = val.(optionstr{cont}{cont2});
        end
        if val
            namefunction = [namefunction '_' optionstr{cont}{1}];
            updatenamefunction = [updatenamefunction '_' optionstr{cont}{1}];
        end
    end
    PGDalgorithm.function.mapping{idim} = str2func([namefunction '_' namedimFile]);
    PGDalgorithm.update.function{idim} = str2func([updatenamefunction '_' namedimFile]);
end

%% NON-CONSTANT COEFFICIENTS 

disp('COMPUTATION PGD COEFFICIENTS')
if loadPGDfromFile(1)
    disp(['PGD data loaded from file ' PGDfile]),
    load([PGDfile '_PGDcoefs.mat'])
else
    dimKTHETA = findPGDdimension('KTHETA',parameters.PGDdimensions);
    if dimKTHETA
        parametersAUX = parameters;
        parametersAUX.PGDdimensions{dimKTHETA,1} = 'K';
        PGDmeshesAUX = PGDmeshes;
        PGDmeshesAUX(dimKTHETA).X = PGDmeshesAUX(dimKTHETA).X(1:parametersAUX.meshes(dimKTHETA).nK+1,1);
        PGDcoefs = computePGDcoefs(parametersAUX,PGDmeshesAUX);
        clear parametersAUX PGDmeshesAUX
    else
        PGDcoefs = computePGDcoefs(parameters,PGDmeshes);
    end
end

%% INCIDENT WAVE

if ~loadPGDfromFile(2)
    disp('COMPUTATION PGD INCIDENT WAVE')
    if ~isempty(parameters.iwparam.load)
        disp(['PGDip loaded from file ' parameters.iwload])
        load(parameters.iwparam.load)
    else
        PGDip = computePGDincidentWave(PGDmeshes,PGDalgorithm,parameters);
        if ~isempty(parameters.iwparam.save)
            disp(['PGDip saved in file ' parameters.iwsave])
            bigsave(parameters.iwparam.save,'PGDip')
        end
    end
end

%% MATRICES AND RHS

disp('COMPUTATION OF PGD MATRICES AND RHS')
if loadPGDfromFile(2)
    disp(['PGD data loaded from file ' PGDfile]),
    load([PGDfile '_PGDmatrices.mat'])
    load([PGDfile '_PGDrhs.mat'])
else
    [PGDmatrices,PGDrhs] = computePGDMatrices(PGDmeshes,PGDcoefs,parameters,PGDalgorithm,PGDip);
    if savePGDintoFile(3) == 2
        PGDmassmatrices = struct();
        for i = 1:parameters.nOfPGDdimensions, PGDmassmatrices(i).M = PGDmatrices(i).M; end
    end
end




