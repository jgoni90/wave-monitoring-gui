function joined_var = MSL_joinVariable(joined_var,var,split_info)

datavar = whos('var');

switch datavar.class
   
    case 'struct'
        error('Join for structures not implemented yet')
        
    case 'cell'
        error('Join for cells not implemented yet')
        
    case 'double'
        
        %Some assignments
        pos     = split_info.file;
        dimstep = split_info.dimstep;
        dimsize = split_info.dimsize;
        ndim    = split_info.ndim;
        dim     = split_info.dim;
        
        %Initialize the joined variable
        if pos == 0
            sizeVar = size(var);
            sizeVar(dim) = dimsize;
            joined_var = zeros(sizeVar);
        end

        %Non split dimensions
        colondims = setdiff(1:ndim,dim);
        S = struct('type','()','subs',{cell(ndim,1)});
        for i = colondims, S.subs{i} = ':'; end
        
        %Join in dimension 'dim'
        S.subs{dim} = pos*dimstep+1:min(dimstep*(pos+1),dimsize);
        
        %Assign joined variable
        joined_var = subsasgn(joined_var,S,var);
        
    otherwise
        error(['Join for variables of class ' datavar.class ' not implemented yet'])
        
end