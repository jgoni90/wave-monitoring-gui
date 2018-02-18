function mysave(file,list_var,varargin)

%Check input
if ischar(list_var), list_var = {list_var}; end
if numel(list_var) ~= numel(varargin)
    error('Size of list_var does not coincide')
end

%Global limit size
global MSL_LIMIT_GB
MSL_LIMIT_GB = 0.5;

%Size in Gb
gb = MSL_getGb(varargin{:});

%Assign the list of variables in the current WS (no copy!)
nOfvariables = numel(varargin);
variables2assign = 1:nOfvariables;
for ivar = variables2assign
    expression = [list_var{ivar} '= varargin{ivar};'];
    eval(expression)
end

%Always split variable from the beginnig (file = 0)
iniFile = 0;
MSL_SAVECOUNTER = iniFile;

%Save variables which dont require splitting
inisave_pos = gb < MSL_LIMIT_GB;
if all(inisave_pos)
    save(file,list_var{inisave_pos});
    return
elseif any(inisave_pos) && iniFile == 0
    run MSL_updateOutputFile
    save(outputfile,list_var{inisave_pos});
else
    %Do nothing, either all the variables should be split or they have been
    %already saved
end

%Split procedure
nOfFiles = ceil(gb / MSL_LIMIT_GB);
cum_nOfFiles_splitVar = cumsum(nOfFiles(~inisave_pos));
[posSplit,posVarFile] = getPositionsInArray(cum_nOfFiles_splitVar,iniFile);
split_var_pos = 1:nOfvariables;
split_var_pos = split_var_pos(~inisave_pos);
for i = split_var_pos(posSplit:end)
    
    inOfFiles = nOfFiles(i);
    for j = posVarFile:inOfFiles-1
        
        %Split
        [split_var,split_info] = MSL_splitVariable(list_var{i},varargin{i},j);
        
        %Save split variable
        run MSL_updateOutputFile
        save(outputfile,'split_var','split_info','MSL_LIMIT_GB');
        clear split_var split_info
    end
    
    posVarFile = 0;
        
end

function [p1,p2] = getPositionsInArray(var,p0)

pos = find(p0-var >= 0,1);
if isempty(pos)
    p1 = 1;
    p2 = max(p0-1,0);
elseif var(pos)-p0 == 0
    p1 = pos;
    p2 = p0-1;
else
    p1 = pos+1;
    p2 = p0 - (var(p1)-var(p1-1)) - 1;
end
    






