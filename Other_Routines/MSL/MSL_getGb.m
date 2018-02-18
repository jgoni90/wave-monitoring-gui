function gb = MSL_getGb(varargin)

%gb = MSL_getGb(var) returns the size "gb" in Gigabytes of list of variables
%given in "varargin"

nOfvariables = numel(varargin);
gb = zeros(nOfvariables,1);

for i = 1:nOfvariables
    ivar = varargin{i};
    varinfo = whos('ivar');
    gb(i) = varinfo.bytes / (1024^3);
end