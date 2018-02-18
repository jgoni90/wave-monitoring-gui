function meshData = MeshCuboidComplete(lims,nela)

% Creates a cuboid structured mesh
% lims  2 * ndim

ndim = size(lims,2);
if numel(nela)==1
    nela = nela*ones(ndim,1);
end
h   = diff(lims)./nela.';
nn  = nela+1;

% elements and corresponding edge/face belonging to each boundary
switch ndim
    case 1
        x   = (lims(1)+(0:nela)*h).';                   % x-coordinates m
        y   = 0*x;                                      % y-coordinates m
        z   = y;                                        % z-coordinates m
        X   = x;                                        % X matrix      m
        
        aux = (1:nela).';
        T   = [aux   aux+1];                            % T matrix
        
        bl = {'b1' 'b2'};     % boundary string labels
        
        bel.(bl{1}) = [1 1];
        bel.(bl{2}) = [size(T,1) 2];
        
        Tb(bl{1})  = 1;
        Tb(bl{2})  = size(X,1);
    case 2
        x      = (lims(1,1)+h(1)*(0:nela(1)).');        % x-coordinates m
        xaux   = x*ones(1,nn(2));
        y      = (lims(1,2)+h(2)*(0:nela(2)).');        % y-coordinates m
        z      = 0*x;                                   % z-coordinates m
        yaux   = ones(nn(1),1)*y.';
        X      = [xaux(:) yaux(:)];                     % X matrix      m
        
        aux    = (1 : (nn(2)-1)*nn(1)).';
        aux(nn(1):nn(1):end) = [];
        T      = [aux aux+1 aux+1+nn(1) aux+nn(1)];     % T matrix
        
        bl = {'b1' 'b2' 'b3' 'b4'};     % boundary string labels
        
        bel.(bl{1}) = [(1:nela(1)).' 1*ones(nela(1),1)];
        bel.(bl{2}) = [nela(1)*(1:nela(2)).' 2*ones(nela(2),1)];
        bel.(bl{3}) = [(size(T,1)-nela(1)+wrev(1:nela(1))).' 3*ones(nela(1),1)];
        bel.(bl{4}) = [wrev(nela(1)*(0:nela(2)-1)).'+1 4*ones(nela(2),1)];
        
        % Connectivity matrix for the boundaries
        nodes = [1 2;2 3;3 4;4 1];
        npe   = size(nodes,2);
        for m=1:length(bl)
            Tb.(bl{m}) = zeros(size(bel.(bl{m}),1),npe);
            for n=1:size(bel.(bl{m}),1)
                Tb.(bl{m})(n,:) = T(bel.(bl{m})(n,1),nodes(bel.(bl{m})(n,2),:));
            end
        end
    case 3
        x    = lims(1,1)+h(1)*(0:nela(1)).';             % x vector      m
        y    = lims(1,2)+h(2)*(0:nela(2)).';             % y vector      m
        z    = lims(1,3)+h(3)*(0:nela(3)).';             % z vector      m
        xaux = x*ones(1,nn(2)*nn(3));
        yaux = ones(nn(1),1)*y.';
        yaux = yaux(:)*ones(1,nn(3));
        zaux = ones(nn(1)*nn(2),1)*z.';
        X    = [xaux(:) yaux(:) zaux(:)];               % X matrix      m
        
        aux = zeros(nn(1)*(nn(2)-1),nn(3)-1);
        for n=1:size(aux,2)
            for m=1:size(aux,1)
                aux(m,n) = (n-1)*nn(1)*nn(2)+m;
            end
        end
        aux(nn(1):nn(1):end,:) = [];
        aux = aux(:);
        T   = [
            aux ...
            aux+1 ...
            aux+1+nn(1) ...
            aux+nn(1) ...
            aux+nn(1)*nn(2) ...
            aux+nn(1)*nn(2)+1 ...
            aux+nn(1)*nn(2)+1+nn(1) ...
            aux+nn(1)*nn(2)+nn(1) ...
            ];
        
        bl  = {'b1' 'b2' 'b3' 'b4' 'b5' 'b6'};
        aux = bsxfun(@plus,(1:nela(1)).',(0:nela(3)-1)*nela(1)*nela(2));
        bel.(bl{1}) = [aux(:) 1*ones(nela(1)*nela(3),1)];
        bel.(bl{3}) = [bel.b1(:,1)+nela(1)*(nela(2)-1) 3*ones(nela(1)*nela(3),1)];
        bel.(bl{4}) = [(1:nela(1):prod(nela)-nela(1)+1).' 4*ones(nela(2)*nela(3),1)];
        bel.(bl{2}) = [bel.b4(:,1)+nela(1)-1 2*ones(nela(2)*nela(3),1)];
        bel.(bl{5}) = [(1:nela(1)*nela(2)).' 5*ones(nela(1)*nela(2),1)];
        bel.(bl{6}) = [bel.b5(:,1)+(nela(3)-1)*nela(1)*nela(2) 6*ones(nela(1)*nela(2),1)];
        
        % Connectivity matrix for the boundaries
        nodes = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;4 3 2 1;5 6 7 8];
        npe   = size(nodes,2);
        for m=1:length(bl)
            belm = bel.(bl{m});
            Tb.(bl{m}) = zeros(size(belm,1),npe);
            for n=1:size(belm,1)
                Tb.(bl{m})(n,:) = T(belm(n,1),nodes(belm(n,2),:));
            end
        end
end
meshData = struct('X',X,'T',T,'x',x,'y',y,'z',z,'nela',nela,'bel',bel,'Tb',Tb);