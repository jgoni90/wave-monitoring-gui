function [res,prev_res1] = computePGDresidual_tensorApproximation...
    (parameters,pgd,meshes,matrices,rhs,algorithm,flag)

%% COMPUTE THE DISCRETE L2 NORM ||rhs - pgd|| / ||pgd||
% NOTE: rhs is a tensor of size d_1 x d_2 x ... x d_n, n = nOfSeparatedDimensions

Uterm = pgd.counters.Uterm;
if isfield(pgd,'alpha'), alpha = pgd.alpha; else alpha = ones(1,Uterm); end

if ~flag
    pos = 1:Uterm;
    if parameters.nOfPGDdimensions == 2
        R = pgd.RB{1}(:,pos)*diag(alpha)*pgd.RB{2}(:,pos).';
    elseif parameters.nOfPGDdimensions == 3
        R = pgd.RB{1}(:,pos)*diag(alpha)*khatrirao(pgd.RB{3}(:,pos),pgd.RB{2}(:,pos)).';
    else
        error('number of separated dimensions not supported in computePGDresidual_tensorApproximation.m')
    end
else
    if parameters.nOfPGDdimensions == 2
        R = pgd.errors.residualMat + pgd.RB{1}(:,Uterm)*pgd.RB{2}(:,Uterm).';
    elseif parameters.nOfPGDdimensions == 3
        R = pgd.errors.residualMat + pgd.RB{1}(:,Uterm)*khatrirao(pgd.RB{3}(:,Uterm),pgd.RB{2}(:,Uterm)).';
    else
        error('number of separated dimensions not supported in computePGDresidual_tensorApproximation.m')
    end
end

res = norm(rhs(:) - R(:)) / norm(R(:));
prev_res1 = R;

