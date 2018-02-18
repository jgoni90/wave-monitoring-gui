function A = tensor_tijkNiNjk(XW,TW,XH,TH,reW,reH,T1)

%Number of elements and number of mesh nodes
nOfElements_W       = size(TW,1);
nOfElements_H       = size(TH,1);
nOfNodes_W          = size(XW,1);
nOfNodes_H          = size(XH,1);
nOfElementNodes_W   = size(TW,2);

%Information of the reference elements
IPw_W    = reW.IPweights1d;
N_W      = reW.N1d;
ngauss_W = length(IPw_W);
IPw_H    = reH.IPweights;
N_H      = reH.N;
Nxi_H    = reH.Nxi;
Neta_H   = reH.Neta;

%Memory allocation
A = zeros(nOfNodes_W,nOfNodes_H);

%Loop in 1D elements
for iElem = 1:nOfElements_W
    Te = TW(iElem,:);
    Xe = XW(Te,:);
    T1e = T1(Te,:);
    
    % Jacobian for straight elements
    Jac = 0.5 * (Xe(2)-Xe(1));
    
    A_e = zeros(nOfElementNodes_W,nOfNodes_H);
    for g = 1:ngauss_W
        N_g = N_W(g,:);
        dline = IPw_W(g)*Jac;
        
        T1e_g = N_g*T1e;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_k = zeros(1,nOfNodes_H);
        for kElem = 1:nOfElements_H
            Te_k = TH(kElem,:);
            Xe_k = XH(Te_k,:);
            
            T1e_gk = T1e_g(Te_k).';
            
            % Jacobian
            Jac_k11 = Nxi_H  *Xe_k(:,1);
            Jac_k21 = Neta_H *Xe_k(:,1);
            Jac_k12 = Nxi_H  *Xe_k(:,2);
            Jac_k22 = Neta_H *Xe_k(:,2);
            detJ = Jac_k11.*Jac_k22 - Jac_k21.*Jac_k12;
            
            % Integration
            f_g = N_H * T1e_gk;
            A_ke = (f_g .* IPw_H .* detJ).' * N_H;
            
            % Assembling
            A_k(:,Te_k) = A_k(:,Te_k) + A_ke;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A_e = A_e + N_g'*A_k*dline;
    end

    % Assembling
    A(Te,:) = A(Te,:) + A_e;
end

