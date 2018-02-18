function [X,T] = MeshCuboid(lims,nela)

% Creates a cuboid structured mesh (p=1)
% lims  2 * ndim: matrix 2 x ndim
% nela must be a column vector

% X is a matrix nNodes * nDim
% T is a matrix nEl    * nNodesPerEl

% By Raul Hospital (2013)
% e-mail: raul.hospital.bravo@gmail.com

ndim = size(lims,2);
if numel(nela)==1
    nela = nela*ones(ndim,1);
end
h   = diff(lims)./nela.';
nn  = nela+1;

switch ndim
    case 1
        X   = (lims(1)+(0:nela)*h).';               % x-coordinates     m
        aux = (1:nela).';
        T   = [aux   aux+1];                        % T matrix
    case 2
        x      = (lims(1,1)+h(1)*(0:nela(1)).');    % x-coordinates     m
        xaux   = x*ones(1,nn(2));
        y      = (lims(1,2)+h(2)*(0:nela(2)).');    % y-coordinates     m
        yaux   = ones(nn(1),1)*y.';
        X      = [xaux(:) yaux(:)];                 % X matrix          m
        aux    = (1 : (nn(2)-1)*nn(1)).';
        aux(nn(1):nn(1):end) = [];
        T      = [aux aux+1 aux+1+nn(1) aux+nn(1)]; % T matrix
    case 3
        x    = lims(1,1)+h(1)*(0:nela(1)).';        % x vector          m
        y    = lims(1,2)+h(2)*(0:nela(2)).';        % y vector          m
        z    = lims(1,3)+h(3)*(0:nela(3)).';        % z vector          m
        xaux = x*ones(1,nn(2)*nn(3));
        yaux = ones(nn(1),1)*y.';
        yaux = yaux(:)*ones(1,nn(3));
        zaux = ones(nn(1)*nn(2),1)*z.';
        X    = [xaux(:) yaux(:) zaux(:)];           % X matrix          m
        aux  = zeros(nn(1)*(nn(2)-1),nn(3)-1);
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
end