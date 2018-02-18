function pgd = initializePGD(algorithm,meshes,matrices,parameters)

%% RECOVER PGD

if parameters.recoverPGD.value
    file = load(parameters.recoverPGD.file);
    filefields = fieldnames(file);
    if any(mystrcmp(filefields,'pgd'))
        pgdfield = 'pgd';
    elseif any(mystrcmp(filefields,'projectedpgd'))
        pgdfield = 'projectedpgd';
    else
        error('recover PGD is not found!!')
    end
    file.(pgdfield).counters.save = 0;
    file.(pgdfield).counters.residual = 0;
    file.(pgdfield).errors.residualMat = 0;
    pgd = file.(pgdfield);
    return
end

%% REDUCED BASIS

pgd = struct();
pgd.RB = cell(parameters.nOfPGDdimensions,1);

nOfTensorNodes = 1;
for idim = 1:parameters.nOfPGDdimensions
    nOfNodesdim = size(meshes(idim).X,1);
    nOfTensorNodes = nOfTensorNodes*nOfNodesdim;
    pgd.RB{idim} = zeros(nOfNodesdim,parameters.maxterms);
end

if algorithm.dual.value
    file = load(algorithm.dual.RBfile);
    pgd.RBdata = file.pgd.RB;
    clear file
end

%% COUNTERS

pgd.counters.Uterm         = 0;
pgd.counters.proj          = 0;
pgd.counters.projTerm      = 1;
pgd.counters.save          = 0;
pgd.counters.residual      = 0;
pgd.counters.residualUterm = 1;

%% UPDATE

if algorithm.update.value
    
    pgd.update = struct();
    updateParameters = algorithm.update.parameter;
    
    for iparam = 1:length(updateParameters)
        param = updateParameters{iparam};
        
        pgd.update.(param).updateI = [];
        pgd.update.(param).updateJ = [];
        pgd.update.(param).updateV = [];
        pgd.update.(param).updatef = [];
    end
    
elseif parameters.residualEachTerm > 0
    
    if strcmpi(parameters.residualType,'DRFS')
        pgd.errors.residualMat = zeros(nOfTensorNodes,1);
    elseif any(strcmpi(parameters.residualType,{'LMC','TENSOR_APP'}))
        pgd.errors.residualMat = 0;
    end
    
end

%% ERRORS

if algorithm.projection.value
    pgd.projection.errors.residual = cell(1,1);
    pgd.errors.residual{1} = zeros(algorithm.projection.projectEachTerm,1);
    pgd.errors.iteration{1} = zeros(algorithm.projection.projectEachTerm,1);
    pgd.errors.nOfiterations{1} = zeros(algorithm.projection.projectEachTerm,1);
else
    pgd.errors.residual{1} = zeros(parameters.maxterms,1);
    pgd.errors.iteration{1} = zeros(parameters.maxterms,1);
    pgd.errors.nOfiterations{1} = zeros(parameters.maxterms,1);
end



