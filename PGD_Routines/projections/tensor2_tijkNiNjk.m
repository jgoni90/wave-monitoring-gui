function A = tensor2_tijkNiNjk(XW,TW,XH,TH,reW,reH,T1)

%Number of tensors
nOfTensors = size(T1,3);

%Number of elements and number of mesh nodes
nOfElements_W       = size(TW,1);
nOfElements_H       = size(TH,1);
nOfNodes_W          = size(XW,1);
nOfNodes_H          = size(XH,1);
nOfElementNodes_W   = size(TW,2);
nOfElementNodes_H   = size(TH,2);

%Information of the reference elements
IPw_W    = reW.IPweights1d;
N_W      = reW.N1d;
ngauss_W = length(IPw_W);
IPw_H    = reH.IPweights;
N_H      = reH.N;
Nxi_H    = reH.Nxi;
Neta_H   = reH.Neta;

%Memory allocation
A = zeros(nOfNodes_W,nOfNodes_H,nOfTensors);

%Loop in 1D elements
for iElem = 1:nOfElements_W
    Te = TW(iElem,:);
    Xe = XW(Te,:);
    T1e = T1(Te,:,:);
    
    % Jacobian for straight elements
    Jac = 0.5 * (Xe(2)-Xe(1));
    
    A_e = zeros(nOfElementNodes_W,nOfNodes_H,nOfTensors);
    for g = 1:ngauss_W
        N_g = N_W(g,:);
        dline = IPw_W(g)*Jac;
        
        T1e_g = zeros(nOfNodes_H,nOfTensors);
        for nt = 1:nOfTensors
            T1e_g(:,nt) = N_g*T1e(:,:,nt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_k = zeros(1,nOfNodes_H,nOfTensors);
        A_ke = zeros(1,nOfElementNodes_H,nOfTensors);
        for kElem = 1:nOfElements_H
            Te_k = TH(kElem,:);
            Xe_k = XH(Te_k,:);
            
            T1e_gk = T1e_g(Te_k,:);
            
            % Jacobian
            Jac_k11 = Nxi_H  *Xe_k(:,1);
            Jac_k21 = Neta_H *Xe_k(:,1);
            Jac_k12 = Nxi_H  *Xe_k(:,2);
            Jac_k22 = Neta_H *Xe_k(:,2);
            detJ = Jac_k11.*Jac_k22 - Jac_k21.*Jac_k12;
            
            % Integration & assembling
            f_g = N_H * T1e_gk;
            for nt = 1:nOfTensors
                A_ke(:,:,nt) = (f_g(:,nt) .* IPw_H .* detJ).' * N_H;
                A_k(:,Te_k,nt) = A_k(:,Te_k,nt) + A_ke(:,:,nt);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for nt = 1:nOfTensors
            A_e(:,:,nt) = A_e(:,:,nt) + N_g'*A_k(:,:,nt)*dline;
        end
    end

    % Assembling
    A(Te,:,:) = A(Te,:,:) + A_e;
end

