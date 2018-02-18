function myload(file,structName)

%Check input
if nargin == 1
    putIntoStruct = false;
elseif nargin == 2
    putIntoStruct = true;
end
    
%Check files
MSL_FILES = dir(['*' file '*MSL*.mat*']);
MSL_ISDIR = [MSL_FILES.isdir]';
if isempty(MSL_FILES)
    var = load(file);
    if putIntoStruct
        assignin('caller',structName,var);
    else
        names = fieldnames(var);
        for i = 1:length(names), assignin('caller',names{i},var.(names{i})); end
    end
    return
elseif all(MSL_ISDIR)
    error([file ' is/are directory/ies, not .mat files'])
elseif any(MSL_ISDIR)
    warning('Some files coincide with directories. The routine may fail')
    MSL_FILES = MSL_FILES(~MSL_ISDIR);
else
    %Do nothing, all files seem to be correct
end
files_names = MSL_sortFiles({MSL_FILES.name}');

%Join procedure
nOfFiles = length(files_names);
joined_var = [];
for i = 1:nOfFiles
    
    var = load(files_names{i});
    names = fieldnames(var);
    if isfield(var,'split_info') %Join the split variable
        
        %Update joined variable (no copy!)
        joined_var = MSL_joinVariable(joined_var,var.split_var,var.split_info);
        
        %Assign when join procedure is finished
        nOfSplitFiles = ceil(var.split_info.dimsize /  var.split_info.dimstep);
        if var.split_info.file == nOfSplitFiles-1
            if putIntoStruct
                oStruct.(var.split_info.name) = joined_var;
            else
                assignin('caller',var.split_info.name,joined_var);
            end
            joined_var = [];
        end
        
    else %Just assign the non-split variables
        if putIntoStruct
            for j = 1:length(names), oStruct.(names{j}) = var.(names{j}); end
        else
            for j = 1:length(names), assignin('caller',names{j},var.(names{j})); end
        end
    end
    
    clear var
        
end

%Output struct if joined variables have been defined
if putIntoStruct
    assignin('caller',structName,oStruct);
end



