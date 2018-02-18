function C = PGDmassMatrix1D(X,T,theReferenceElement,coefficients,var,elemCoef)

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);
nOfBoundaryNodes = length(unique(T));
nOfLinearNodes = nOfElements + 1;
nOfElementNodes = theReferenceElement.degree + 1;
nOfInnerNodes = nOfBoundaryNodes - nOfLinearNodes;

%General memory allocation
NNZ = nOfElementNodes*(2*nOfLinearNodes + nOfInnerNodes);
aux_ones = ones(1,nOfElementNodes);

%Loop in number of coefficients
nOfCoefs = numel(coefficients);
C = cell(nOfCoefs,1);
for icoef = 1:nOfCoefs
    
    %Initialization and memory allocation
    coef = coefficients{icoef}{var};
    nOfTerms = size(coef,2);
    C{icoef} = cell(nOfTerms,1);
    I = zeros(NNZ,1);
    J = I;
    Creal = zeros(NNZ,nOfTerms);
    Cimag = Creal;
    
    %Loop in 1D boundary elements
    for ielem = 1:nOfElements
        Te = T(ielem,:);
        Xe = X(Te,:);
        alpha = elemCoef(ielem);
        coefe = coef(Te,:);
        
        %Elemental matrices
        Ce = elementalMatrix(Xe,theReferenceElement,coefe);
        Ce = alpha*Ce; %alpha is constant element by element
        
        %Assembling
        Te_transp = transpose(Te);
        aux_row = Te_transp(:,aux_ones);
        aux_col = Te(aux_ones,:);
        indK = (ielem-1)*nOfElementNodes^2+1:ielem*nOfElementNodes^2;
        I(indK) = aux_row(:);
        J(indK) = aux_col(:);
        for i = 1:nOfTerms
            iCe = Ce(:,:,i);
            Creal(indK,i) = real(iCe(:));
            Cimag(indK,i) = imag(iCe(:));
        end
    end
    
    %Build sparse matrix
    pos = I > 0;
    for i = 1:nOfTerms
        C{icoef}{i} = sparse(I(pos),J(pos),...
            Creal(pos,i) + sqrt(-1)*Cimag(pos,i),...
            nOfNodes,nOfNodes);
    end
    
end

%______________________________________________________________
function Ce = elementalMatrix(Xe,theReferenceElement,coefe)

nOfTerms = size(coefe,2);
nOfElementNodes = size(Xe,1);
Ce = zeros(nOfElementNodes,nOfElementNodes,nOfTerms);

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
    %Contribution of the current integration point to the elemental matrix
    elemN = (N_g')*N_g*dline;
    for i = 1:nOfTerms, Ce(:,:,i) = Ce(:,:,i) + coef(i)*elemN; end
end


