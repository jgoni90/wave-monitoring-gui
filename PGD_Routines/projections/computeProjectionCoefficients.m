function [alpha,Mat,f] = computeProjectionCoefficients...
            (pgd,mode,prevMat,prevf,matrices,parameters,pgdRB)

%Initialize
Uterm                   = pgd.counters.Uterm;        
Mat                     = ones(Uterm+1,Uterm+1); %diag(Mat) = (mode,mode)_L2 = 1
f                       = ones(Uterm+1,1);
f(1:Uterm)              = prevf;
Mat(1:Uterm,1:Uterm)    = prevMat;
vstored                 = cell(1,parameters.nOfPGDdimensions);
for j = 1:parameters.nOfPGDdimensions, vstored{j} = matrices(j).M * mode{j}; end

%Compute last column -but diagonal term- of matrix Mat as (mode,pgd.RB(:)(:,1:Uterm))_L2
for i = 1:Uterm
    for j = 1:parameters.nOfPGDdimensions
        Mat(i,Uterm+1) = Mat(i,Uterm+1) * (pgd.RB{j}(:,i)' * vstored{j});
    end
end

%Complete the hermitian system matrix Mat
Mat(Uterm+1,1:Uterm) = Mat(1:Uterm,Uterm+1)';

%Last term of separable RHS vector computed with pgdRB structure
v = 1;
for j = 1:parameters.nOfPGDdimensions
    v = v .* ((mode{j}'*matrices(j).M) * pgdRB{j});
end
f(Uterm+1) = sum(v);

%Compute alpha coefficients
alpha = (Mat \ f).';