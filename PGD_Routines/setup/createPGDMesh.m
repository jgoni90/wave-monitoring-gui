function [X,T] = createPGDMesh(meshType,meshParameters)

switch meshType
    case {'K','THETA'}
        
        %Parameters
        nDeg = meshParameters.nDeg;
        nOfNodes = meshParameters.ininOfNodes;
        if strcmpi(meshType,'K')
            a = 2*pi/meshParameters.Tend;
            b = 2*pi/meshParameters.Tini;
        elseif strcmpi(meshType,'THETA')
            a = meshParameters.THETAini;
            b = meshParameters.THETAend;
        end

        %Nodal positions
        nElems = round((nOfNodes-1)/nDeg); 
        nOfNodes = nElems*nDeg + 1; 
        X = chev(nOfNodes,a,b);
        
        %Connectivity
        T = create1Dconec(nElems,nDeg); 
        
    case 'KTHETA'
        
        %Parameters
        fieldNames = fieldnames(meshParameters);
        nDeg = meshParameters.nDeg;
        if nDeg > 1, error('KTHETA meshes for degree > 1 not implemented yet'), end
        if all(mystrcmp({'Tend','Tini'},fieldNames))
            wini = 2*pi/meshParameters.Tend;
            wend = 2*pi/meshParameters.Tini;
        elseif all(mystrcmp({'Kini','Kend'},fieldNames))
            wini = meshParameters.Kini;
            wend = meshParameters.Kend;
        end
        THETAini = meshParameters.THETAini;
        THETAend = meshParameters.THETAend;
        nTHETA = meshParameters.nTHETA;
        equalSpaced = meshParameters.equalSpaced;
        nw = meshParameters.nK;
        
        %2D Mesh
        [Xaux,T] = MeshCuboid([wini,THETAini ; wend,THETAend],[nw ; nTHETA]);
        
        %Nodal positions
        if equalSpaced
            X = Xaux;
        else
            XW = chev(nw+1,wini,wend);
            XTHETA = chev(nTHETA+1,THETAini,THETAend);
            X = mesh1Dto2D(XW,XTHETA);
        end
        
    case 'DoQUAD'

        %Parameters
        x = meshParameters.x;
        x1 = x(1); x2 = x(end);
        y = meshParameters.y;
        y1 = x(1); y2 = x(end);
        nx = numel(x) - 1;
        ny = numel(y) - 1;

        %2D Mesh
        [~,T] = MeshCuboid([x1,y1 ; x2,y2],[nx ; ny]);
        X = mesh1Dto2D(x,y);
end


%% CHEBYSHEV DISTRIBUTION

function X = chev(n,a,b)

X = zeros(n,1);
m = n - 2;
for i = 1:m
    X(n-i) = 0.5*(a + b) + 0.5*(b - a)*cos((2*i-1)*pi/(2*m));
end
X(1) = a;
X(n) = b;


%% 2D MESH FROM 1D MESHES

function X = mesh1Dto2D(x,y)

nx = size(x,1); 
ny = size(y,1);
ini = 1;
fin = nx;
X = zeros(nx*ny,2);
for i = 1:ny
    X(ini:fin,:) = [x y(i)*ones(size(x))];
    ini = fin + 1;
    fin = fin + nx;
end



