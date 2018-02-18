function f = projectionReUIVectorX_NonSeparable_K(meshes,mode)

% X1,T1,re1:        K mesh
% X2,T2,re2:        X mesh
% X3,T3,re3:        THETA mesh
% R(X),T(THETA):    nodal vectors for X and THETA dimensions to be integrated

X1  = meshes(2).X;
X2  = meshes(1).X;
X3  = meshes(3).X;
T1  = meshes(2).T;
T2  = meshes(1).T;
T3  = meshes(3).T;
re1 = meshes(2).referenceElement;
re2 = meshes(1).referenceElement;
re3 = meshes(3).referenceElement;
R   = mode{1};
T   = mode{2};

%Vectorized elements
nev1 = 100;
nev2 = 100;
nev3 = 100;

%Number of elements and number of mesh nodes
[ne1,nen1] = size(T1);
[ne2,nen2] = size(T2);
[ne3,nen3] = size(T3);

%Information of the reference element
IPw1 = re1.IPweights1d; ng1 = length(IPw1);
IPw2 = re2.IPweights1d; ng2 = length(IPw2);
IPw3 = re3.IPweights1d; ng3 = length(IPw3);
N1 = re1.N1d;
N2 = re2.N1d;
N3 = re3.N1d;

%Precomputed tensors in the isoparametric space
N_i_g1 = bsxfun(@times,N1',reshape(IPw1,1,ng1));
N_g2_l = bsxfun(@times,N2,reshape(IPw2,ng2,1));
N_g3_m = bsxfun(@times,N3,reshape(IPw3,ng3,1));
N_g1g2g3_ilm = bsxfun(@times,reshape(N1,ng1,1,1,nen1,1,1),...
    reshape(N2,1,ng2,1,1,nen2,1));
N_g1g2g3_ilm = bsxfun(@times,N_g1g2g3_ilm,reshape(N3,1,1,ng3,1,1,nen3));
N_g1g2g3_ilm = reshape(N_g1g2g3_ilm,ng1*ng2*ng3,nen1*nen2*nen3);

%Memory allocation
f = zeros(nen1*ne1,1);
IF = f;

%Vectors of vectorized elements
ve1 = 1:nev1:ne1; if ve1(end) ~= ne1, ve1 = [ve1,ne1]; end, ve1 = ve1(2:end);
ve2 = 1:nev2:ne2; if ve2(end) ~= ne2, ve2 = [ve2,ne2]; end, ve2 = ve2(2:end);
ve3 = 1:nev3:ne3; if ve3(end) ~= ne3, ve3 = [ve3,ne3]; end, ve3 = ve3(2:end);

%Useful transformations
T1_t = T1';
X1_t = X1';
T2_t = T2';
X2_t = X2';
T3_t = T3';
X3_t = X3';

%Loop in vectorized elements
inie1 = 1;
for ie1 = ve1
    
    %Coordinates of elements
    ivE1    = inie1:ie1;
    nE1     = length(ivE1);
    iTe1    = T1_t(:,ivE1);                                           
    Te1     = reshape(iTe1,1,nE1*nen1);                     
    K       = reshape(X1_t(:,Te1),nen1,1,1,nE1,1,1);
    
    %Jacobian (straight elements!!)
    J1 = (1/2) * abs((K(end,:,:,:,:,:) - K(1,:,:,:,:,:)));
    J1 = reshape(J1,1,1,1,nE1,1,1);
    
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
        X       = reshape(X2_t(:,Te2),1,nen2,1,1,nE2,1);
        
        %Jacobian (straight elements!!)
        J2 = (1/2) * abs((X(:,end,:,:,:,:) - X(:,1,:,:,:,:)));
        
        %Nodal coordinates of function S and product with jacobian
        iR = bsxfun(@times,reshape(R(Te2),1,nen2,1,1,nE2,1),J2);
        iR = N_g2_l * reshape(iR,nen2,nE2);
        
        inie3 = 1;
        for ie3 = ve3
            
            %Coordinates of elements
            ivE3    = inie3:ie3;
            nE3     = length(ivE3);
            iTe3    = T3_t(:,ivE3);                                           
            Te3     = reshape(iTe3,1,nE3*nen3);                     
            THETA   = reshape(X3_t(:,Te3),1,1,nen3,1,1,nE3);
            
            %Jacobian (straight elements!!)
            J3 = (1/2) * abs((THETA(:,:,end,:,:,:) - THETA(:,:,1,:,:,:)));
            
            %Nodal coordinates of function T and product with jacobian
            iT = bsxfun(@times,reshape(T(Te3),1,1,nen3,1,1,nE3),J3);
            iT = N_g3_m * reshape(iT,nen3,nE3);
            
            %Incident wave valued at coordinates of the nodes
            UI = bsxfun(@times,K,cos(THETA));
            UI = bsxfun(@times,X,UI);
            UI = cos(UI);
            
            %Interpolation of incident wave and product with jacobian
            fval = N_g1g2g3_ilm * reshape(UI,nen1*nen2*nen3,nE1*nE2*nE3);
            fval = reshape(fval,ng1,ng2,ng3,nE1,nE2,nE3);
            fval = bsxfun(@times,fval,J1);
            
            %Vectorized elemental contribution to the vector
            tensor_beta = bsxfun(@times,fval,reshape(iT,1,1,ng3,1,1,nE3));
            tensor_beta = bsxfun(@times,tensor_beta,reshape(iR,1,ng2,1,1,nE2,1));
            tensor_beta = sum(tensor_beta,2); %contraction in ng2
            tensor_beta = sum(tensor_beta,3); %contraction in ng3
            tensor_beta = sum(tensor_beta,5); %contraction in nen2
            tensor_beta = sum(tensor_beta,6); %contraction in nen3
            fe = N_i_g1 * reshape(tensor_beta,ng1,nE1);
            
            % Assembling
            f(indF) = f(indF) + fe(:);
            
            inie3 = ie3 + 1;
        end
    
        inie2 = ie2 + 1;
    end
    
    inie1 = ie1 + 1;
end

%Output
f = accumarray(IF,f);
    
    












