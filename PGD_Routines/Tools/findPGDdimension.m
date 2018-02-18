function pos = findPGDdimension(param,dimensions)

nOfDimensions = size(dimensions,1);
    
if strcmp(param,'PARAMETRICDIM')
    pos = [];
    for i = 1:nOfDimensions
        idim = dimensions(i,:);
        if idim{4}
            pos = [pos i];
        end
    end
    
elseif strcmp(param(end),'*')
    pos = [];
    paramcopy = param(1:end-1);
    nCharparam = length(paramcopy);
    for i = 1:nOfDimensions
        idim = dimensions(i,:);
        iname = idim{1};
        try if strcmp(iname(1:nCharparam),paramcopy), pos = [pos i]; end
        catch, continue, end
    end

else
    if ~iscell(param), paramcopy = {param}; else paramcopy = param; end
    pos = zeros(size(paramcopy));
    for i = 1:nOfDimensions
        idim = dimensions(i,:);
        cond = strcmpi(idim{1},paramcopy);
        if any(cond)
            pos(cond) = i;
        end
    end
    
end

if isempty(pos), pos = 0; end
        