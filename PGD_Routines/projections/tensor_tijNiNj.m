function A = tensor_tijNiNj(XW,TW,XH,TH,reW,reH,T1)

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
Nxi_W    = reW.N1dxi;
ngauss_W = length(IPw_W);
IPw_H    = reH.IPweights1d;
N_H      = reH.N1d;
Nxi_H    = reH.N1dxi;
ngauss_H = length(IPw_H);

%Memory allocation
A = zeros(nOfNodes_W,nOfNodes_H);

%Loop in 1D elements
for iElem = 1:nOfElements_W
    Te = TW(iElem,:);
    Xe = XW(Te,:);
    T1e = T1(Te,:);

    A_e = zeros(nOfElementNodes_W,nOfNodes_H);
    for g = 1:ngauss_W
        N_g = N_W(g,:);
        Nxi_g = Nxi_W(g,:);
        xyDer_g = Nxi_g*Xe;
        xyDerNorm_g = norm(xyDer_g);
        dline = IPw_W(g)*xyDerNorm_g;
        
        T1e_g = N_g*T1e;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_k = zeros(1,nOfNodes_H);
        for kElem = 1:nOfElements_H
            Te_k = TH(kElem,:);
            Xe_k = XH(Te_k);
            
            T1e_gk = T1e_g(Te_k).';
            
            A_ke = zeros(1,nOfElementNodes_H);
            for gk = 1:ngauss_H
                N_gk = N_H(gk,:);
                Nxi_gk = Nxi_H(gk,:);
                xyDer_gk = Nxi_gk*Xe_k;
                xyDerNorm_gk = norm(xyDer_gk);
                dline_k = IPw_H(gk)*xyDerNorm_gk;
                
                f_g = N_gk*T1e_gk;
                
                A_ke = A_ke + f_g*N_gk*dline_k;
            end
            
            % Assembling
            A_k(:,Te_k) = A_k(:,Te_k) + A_ke;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A_e = A_e + N_g'*A_k*dline;
    end

    % Assembling
    A(Te,:) = A(Te,:) + A_e;
end

