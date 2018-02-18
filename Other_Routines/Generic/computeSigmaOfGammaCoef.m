function f = computeSigmaOfGammaCoef(X,T,sigma)

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Memory allocation
f = zeros(nOfNodes,1);
gamma = zeros(nOfNodes,1);

%Loop in 1D boundary elements
for ielem = 1:nOfElements
    Te = T(ielem,:);
    Xe = X(Te,:);
    sigmae = sigma(Te,:);
    INDEX_PML_BOUNDARY = getIndex(Xe);
    gamma(Te) = sigmae(:,INDEX_PML_BOUNDARY);
end

function i = getIndex(X)

tol = 1e-5;
j = 1;
res = tol/2;
while res < tol && j < size(X,1)
    res = abs((X(j,1) - X(j+1,1)) / X(j+1,1));
    j = j+1;
end
if res < tol, i = 2; else i = 1; end
