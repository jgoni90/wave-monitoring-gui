function [err,prev_err] = computePGDerror(varargin)

%PGD parameters should be the first input

parameters = varargin{1};
switch parameters.residualType
    case {'LMC','lmc'}
        [err,prev_err] = computePGDresidualContribution(varargin{:});
    case {'DRFS','drfs'}
        [err,prev_err] = computePGDresidual(varargin{:});
    case {'TENSOR_APP','tensor_app'}
        [err,prev_err] = computePGDresidual_tensorApproximation(varargin{:});
    otherwise
        err = -1;
        prev_err = -1;
end

