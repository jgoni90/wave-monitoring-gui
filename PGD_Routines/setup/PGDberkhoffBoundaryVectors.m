function [fx,fy] = PGDberkhoffBoundaryVectors(X,T2D,T,theReferenceElement,elements2D,derN2Don1D,u,coef1,coef2)

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Memory allocation
NNZ = numel(unique(T));
nOfTermsu = size(u,2);
nOfTerms1 = size(coef1,2);
nOfTerms2 = size(coef2,2);
fx = cell(nOfTerms1,1);
fy = cell(nOfTerms2,1);
fx(:) = {spalloc(nOfNodes,nOfTermsu,NNZ)};
fy(:) = {spalloc(nOfNodes,nOfTermsu,NNZ)};
 
%Loop in 1D boundary elements
for ielem = 1:nOfElements
    Te = T(ielem,:);
    Xe = X(Te,:);
    derN2De = derN2Don1D(:,:,:,ielem);
    coef1e = coef1(Te,:);
    coef2e = coef2(Te,:);
    elem2D = elements2D(ielem);
    ue = u(T2D(elem2D,:),:);
    
    %Elemental matrices (Sparse assembling to avoid heavy allocation... can be slow)
    for i = 1:nOfTerms1
        fe = elementalVector(Xe,theReferenceElement,coef1e(:,i),ue,derN2De(:,:,1),1);
        fx{i}(Te,:) = fx{i}(Te,:) + fe;
    end
    for i = 1:nOfTerms2
        fe = elementalVector(Xe,theReferenceElement,coef2e(:,i),ue,derN2De(:,:,2),2);
        fy{i}(Te,:) = fy{i}(Te,:) + fe;
    end
end


%______________________________________________________________
function fe = elementalVector(Xe,theReferenceElement,coefe,ue,derN,var)

nOfTerms = size(ue,2);
nOfElementNodes = size(Xe,1);
fe = zeros(nOfElementNodes,nOfTerms);

%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

%Number of Gauss points
ngauss = length(IPw);

%Loop in integration points
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    %Integration weight
    xyDer_g = Nxi_g*Xe;
    xyDerNorm_g = norm(xyDer_g);
    dline=IPw(g)*xyDerNorm_g;
    %Variable coefficients (vectorized)
    coef = N_g*coefe;
    uder = derN(g,:)*ue;
    %Unit normal to the boundary
    t_g = xyDer_g/xyDerNorm_g;
    n_g = [t_g(2) -t_g(1)];
    %Contribution of the current integration point to the elemental matrix
    elemN = n_g(var)*coef*(N_g')*dline;
    for i = 1:nOfTerms, fe(:,i) = fe(:,i) + uder(i)*elemN; end
end


