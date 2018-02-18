function [res,prev_res1] = computePGDresidualContribution(parameters,pgd,meshes,matrices,rhs,algorithm,flag)

% For multidimensional problem with D dimensions (x_1,...,x_D):
% Compute e{n} = ||Phi{n} / u{n-1}||_L2 for u{n} = sum_{i=1}^n Phi{i}
% with Phi{i} = prod_{k=1}^D F_i^k(x_D)

%% COMPUTE THE CONTRIBUTION OF THE LAST MODE OF THE SEPARATED REPRESENTATION

Uterm = pgd.counters.Uterm;
if Uterm == 1
    res = 1;
    prev_res1 = 0;
    return
end

if isfield(pgd,'alpha'), alpha = pgd.alpha; else alpha = ones(1,Uterm); end

if ~flag %Compute the complete series ||u{n-1}||
    
    res1 = 0;
    if algorithm.projection.projCoeff
        for i = 1:Uterm-1
            for j = 1:Uterm-1
                res1 = res1 + alpha(i) * alpha(j)';
            end
        end
    else
        for i = 1:Uterm-1
            for j = 1:Uterm-1
                prod = alpha(i) * alpha(j)';
                for k = 1:parameters.nOfPGDdimensions
                    prod = prod * pgd.RB{k}(:,j)' * matrices(k).M * pgd.RB{k}(:,i);
                end
                res1 = res1 + prod;
            end
        end
    end
    
else % Compute the series ||u{n-1}|| with recursive formulae:
     % ||u{n-1}|| = ||u{n-2}|| + p_{n-1,n-1} + sum_{i=1}^{n-2} (p_{i,n-1} + conj(p_{i,n-1}))
     % where p_{i,j} = prod_{k=1}^D (F_j^k)' * M * F_i^k, with M = mass matrix 
    
    prev_res1 = pgd.errors.residualMat;
     
    S_i_n2 = 0;
    for i = 1:Uterm-2
        p_i_n2 = 1;
        for j = 1:parameters.nOfPGDdimensions
            p_i_n2 = p_i_n2 * pgd.RB{j}(:,Uterm-1)' * matrices(j).M * pgd.RB{j}(:,i);
        end
        S_i_n2 = S_i_n2 + p_i_n2 + p_i_n2';
    end

    p_n1_n1 = 1;
    for j = 1:parameters.nOfPGDdimensions
        p_n1_n1 = p_n1_n1 * pgd.RB{j}(:,Uterm-1)' * matrices(j).M * pgd.RB{j}(:,Uterm-1);
    end

    res1 = prev_res1 + p_n1_n1 + S_i_n2;
end

%Compute ||Phi{n}||
p_n_n = alpha(Uterm) * alpha(Uterm)';
for j = 1:parameters.nOfPGDdimensions
    p_n_n = p_n_n * pgd.RB{j}(:,Uterm)' * matrices(j).M * pgd.RB{j}(:,Uterm);
end

%Compute e{n} and output ||u{n-1}|| for next recursive step if needed
res = real(p_n_n) / real(res1);
prev_res1 = res1;


