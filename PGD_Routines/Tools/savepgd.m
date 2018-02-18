function savepgd(filename,mode,listOfvariables,keyfields,varargin)

% This version DOES NOT ALLOW to store multidimensional structures and
% non-vector cells

POS = 1:mode;
struc2save = getVariables(POS,listOfvariables,keyfields,varargin);
disp('Saving...')
bigsave(filename,'-struct','struc2save')

function struc2save = getVariables(POS,listOfvariables,keyfields,varargin)

for i = 1:size(listOfvariables,1)
    
    if ~listOfvariables{i,2} || any(strcmp(listOfvariables{i,1},keyfields)) || isempty(varargin{:}{i})
        struc2save.(listOfvariables{i,1}) = varargin{:}{i}; %save complete variable
        
    elseif isstruct(varargin{:}{i})
        fields = fieldnames(varargin{:}{i});
        for j = 1:length(fields)
            jlist = {fields{j} 1};
            auxstruct = getVariables(POS,jlist,keyfields,{varargin{:}{i}.(fields{j})});
            struc2save.(listOfvariables{i,1}).(fields{j}) = auxstruct.(fields{j});
        end
        
    elseif iscell(varargin{:}{i})
        for j = 1:length(varargin{:}{i})
            auxstruct = getVariables(POS,listOfvariables(i,:),keyfields,varargin{:}{i}(j));
            struc2save.(listOfvariables{i,1}){j} = auxstruct.(listOfvariables{i,1});
        end
        
    elseif isscalar(varargin{:}{i})
        struc2save.(listOfvariables{i,1}) = varargin{:}{i};
        
    elseif isvector(varargin{:}{i})
        struc2save.(listOfvariables{i,1}) = varargin{:}{i}(POS);
        
    else
        struc2save.(listOfvariables{i,1}) = varargin{:}{i}(:,POS);
    end
    
end

        
    
    