function u = computeFEMatSamplePoints(dataGUI,interpX,parameters)

load(dataGUI)
addpath(genpath([pwd '/../../FEM']))
addpath(genpath([pwd '/../../Callbacks']))

param = findPGDdimension('PARAMETRICDIM',parameters.PGDdimensions);
paramnames = parameters.PGDdimensions(param,1);
k_index = findPGDdimension('K',parameters.PGDdimensions);
theta_index = findPGDdimension('THETA',parameters.PGDdimensions);

nOfParams = numel(paramnames);
nOfSamples = 1; for i = 1:nOfParams, nOfSamples = nOfSamples*size(interpX{param(i)},1); end
u = zeros(size(data.mesh.X,1),nOfSamples);

firstParamSize = size(interpX{param(1)},1);

i = 1;
j = 1;
for isample = 1:nOfSamples
    
    %Parametric data
    if all(mystrcmp(paramnames,{'K' 'THETA'},true))
        k = interpX{k_index}(i);
        bottom = parameters.meshes(k_index).bottomValue;
        theta = interpX{theta_index}(j);
    elseif mystrcmp(paramnames,'KTHETA',true)
        error('KTHETA dim not implemented yet in computeFEMatSamplePoints')
    elseif mystrcmp(paramnames,'THETA',true)
        k = parameters.fixedParametricDims(1);
        bottom = 1;
        theta = interpX{theta_index}(i);
    elseif mystrcmp(paramnames,'K',true)
        k = interpX{k_index}(i);
        bottom = parameters.meshes(k_index).bottomValue;
        theta = parameters.fixedParametricDims(2);
    end
    data.ip.period = 2*pi/sqrt(k*9.81*tanh(k*bottom));
    if isfield(data.bottom,'ccg'), data.bottom = rmfield(data.bottom,'ccg'); end
    data.ip.direction = theta*180/pi;
    
    %FEM solution
    data = run_domain(data,[]);
    data = run_boundary(data,[]);
    u(:,isample) = data.solution;
    
    %Update counters
    if nOfParams == 2
        if j == firstParamSize, i = i + 1; j = 1; else j = j + 1; end
    else
        i = i + 1;
    end
end

