function [split_var,split_info] = MSL_splitVariable(varname,var,pos)

datavar = whos('var');

switch datavar.class
   
    case 'struct'
        error('Splitting for structures not implemented yet')
        
    case 'cell'
        error('Splitting for cells not implemented yet')
        
    case 'double'
        
        %Get dimensions
        if datavar.complex, factor = 2; else factor = 1; end
        [dim,ndim,dimstep,dimsize] = getDim(var,factor);
        split_info = struct('name',varname,'dim',dim,'ndim',ndim,...
            'dimstep',dimstep,'dimsize',dimsize,'file',pos);
        
        %Non splitt dimensions
        colondims = setdiff(1:ndim,dim);
        S = struct('type','()','subs',{cell(ndim,1)});
        for i = colondims, S.subs{i} = ':'; end
        
        %Splitting in dimension 'dim'
        S.subs{dim} = pos*dimstep+1:min(dimstep*(pos+1),dimsize);
        
        %Assign split variable
        split_var = subsref(var,S);
        
    otherwise
        error(['Splitting for variables of class ' datavar.class ' not implemented yet'])
        
end

function [d,nd,s,ds] = getDim(var,f)

global MSL_LIMIT_GB

sizevar = size(var);
nd = numel(sizevar);
[~,d] = max(sizevar);
ds = sizevar(d);
sizevar(d) = [];
s = floor(MSL_LIMIT_GB*(1024^3) / (f*8*prod(sizevar)));