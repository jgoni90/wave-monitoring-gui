function f = nonlin1(meshes,mode)

% X1,T1,re1:        THETA mesh
% X2,T2,re2:        C mesh
% S(C)              nodal vectors for C dimension to be integrated
% This rotine can be used also to build the vector for the C mesh

X1  = meshes(1).X;
X2  = meshes(2).X;
T1  = meshes(1).T;
T2  = meshes(2).T;
re1 = meshes(1).referenceElement;
re2 = meshes(2).referenceElement;
S   = mode{1};

%Vectorized elements
nev1 = 100;
nev2 = 100;

%Number of elements and number of mesh nodes
[ne1,nen1] = size(T1);
[ne2,nen2] = size(T2);

%Information of the reference element
IPw1 = re1.IPweights1d; ng1 = length(IPw1);
IPw2 = re2.IPweights1d; ng2 = length(IPw2);
N1 = re1.N1d;
N2 = re2.N1d;

%Precomputed tensors in the isoparametric space
N_i_g1 = bsxfun(@times,N1',reshape(IPw1,1,ng1));
N_g2_l = bsxfun(@times,N2,reshape(IPw2,ng2,1));
N_g1g2_il = bsxfun(@times,reshape(N1,ng1,1,nen1,1),...
    reshape(N2,1,ng2,1,nen2));
N_g1g2_il = reshape(N_g1g2_il,ng1*ng2,nen1*nen2);

%Memory allocation
f = zeros(nen1*ne1,1);
IF = f;

%Vectors of vectorized elements
ve1 = 1:nev1:ne1; if ve1(end) ~= ne1, ve1 = [ve1,ne1]; end, ve1 = ve1(2:end);
ve2 = 1:nev2:ne2; if ve2(end) ~= ne2, ve2 = [ve2,ne2]; end, ve2 = ve2(2:end);

%Useful transformations
T1_t = T1';
X1_t = X1';
T2_t = T2';
X2_t = X2';

%Loop in vectorized elements
inie1 = 1;
for ie1 = ve1
    
    %Coordinates of elements
    ivE1    = inie1:ie1;
    nE1     = length(ivE1);
    iTe1    = T1_t(:,ivE1);                                           
    Te1     = reshape(iTe1,1,nE1*nen1);                     
    X       = reshape(X1_t(:,Te1),nen1,1,nE1,1);
    
    %Jacobian (straight elements!!)
    J1 = (1/2) * abs((X(end,:,:,:) - X(1,:,:,:)));
    J1 = reshape(J1,1,1,nE1,1);
    
    %Assembling variables
    indF = (inie1-1)*nen1+1:ie1*nen1;
    IF(indF) = Te1;
    
    inie2 = 1;
    for ie2 = ve2
        
        %Coordinates of elements
        ivE2    = inie2:ie2;
        nE2     = length(ivE2);
        iTe2    = T2_t(:,ivE2);                                           
        Te2     = reshape(iTe2,1,nE2*nen2);                     
        K       = reshape(X2_t(:,Te2),1,nen2,1,nE2);
        
        %Jacobian (straight elements!!)
        J2 = (1/2) * abs((K(:,end,:,:) - K(:,1,:,:)));
        
        %Nodal coordinates of function S and product with jacobian
        iS = bsxfun(@times,reshape(S(Te2),1,nen2,1,nE2),J2);
        iS = N_g2_l * reshape(iS,nen2,nE2);
        
        %Non-liner function valued at coordinates of the nodes
        COEF = bsxfun(@times,X,K);
        COEF = cos(COEF);

        %Interpolation of the non-linear function and product with jacobian
        fval = N_g1g2_il * reshape(COEF,nen1*nen2,nE1*nE2);
        fval = reshape(fval,ng1,ng2,nE1,nE2);
        fval = bsxfun(@times,fval,J1);
        
        %Vectorized elemental contribution to the vector
        tensor_beta = bsxfun(@times,fval,reshape(iS,1,ng2,1,nE2));
        tensor_beta = sum(tensor_beta,2); %contraction in ng2
        tensor_beta = sum(tensor_beta,4); %contraction in nen2
        fe = N_i_g1 * reshape(tensor_beta,ng1,nE1);

        % Assembling
        f(indF) = f(indF) + fe(:);
    
        inie2 = ie2 + 1;
    end
    
    inie1 = ie1 + 1;
end

%Output
f = accumarray(IF,f);
    
    












