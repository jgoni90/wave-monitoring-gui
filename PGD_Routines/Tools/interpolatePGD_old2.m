function [u,interpolatedCoef] = interpolatePGD(dim,varargin)

% Return u(x,y) = u(x,y,param1,param2), interpolated at parameter values param1 
% and param2, from reduced basis alpha*R(x,y)*S(param1)*T(param2).
%
%       u = interpolatePGD(dim=4,param1,param2,S,T,R,X_s,T_s,nDeg_s,X_t,T_t,nDeg_t,alpha)
%       u = interpolatePGD(dim=3,param1,S,R,X_s,T_s,nDeg_s,alpha)
%
%   Inputs:
%       dim: number of total dimensions (space + parameters).
%       param: parameter value where the PGD approximation will be
%              interpolated.
%       S,T,R: reduced basis (R(x,y) has to be the last one).
%       X,T: one-dimensional parametric mesh, X: nodal position, T: connectivity.
%       nDeg: degree of the parametric mesh.
%       alpha: column vector of coefficients of the reduced basis.

%Check dimension type
switch dim
    case 3
        pt{1} = varargin{1};
        reducedBasis(1) = varargin(2);
        reducedBasis(2) = varargin(3);
        parametricMesh(1,1) = varargin(4);
        parametricMesh(1,2) = varargin(5);
        nDegParametric = varargin{6};
        nOfParametricDimensions = 1;
    case 4
        pt{1} = varargin{1};
        pt{2} = varargin{2};
        reducedBasis(1) = varargin(3);
        reducedBasis(2) = varargin(4);
        reducedBasis(3) = varargin(5);
        parametricMesh(1,1) = varargin(6);
        parametricMesh(1,2) = varargin(7);
        nDegParametric(1) = varargin{8};
        parametricMesh(2,1) = varargin(9);
        parametricMesh(2,2) = varargin(10);
        nDegParametric(2) = varargin{11};
        nOfParametricDimensions = 2;
end
alpha = varargin{end};

%Look for the non-spacial elements
elem = cell(nOfParametricDimensions,1);
for i = 1:nOfParametricDimensions
    N = nodalConnectivity(parametricMesh{i,2});
    nn = dsearchn(parametricMesh{i,1},pt{i});
    enn = N(nn,:);
    elem{i} = enn(:,1);
    pos = pt{i} - parametricMesh{i,1}(nn) > 0 & enn(:,2) > 0;
    elem{i}(pos) = enn(pos,2);
end
            
%Interpolate non-spacial dimensions (only 1 point)
nOfTerms = size(reducedBasis{1},2);
interpolatedParam = zeros(1,nOfTerms);
interpolatedCoef = alpha;
for i = 1:nOfParametricDimensions
    coord = parametricMesh{i,1}(parametricMesh{i,2}(elem{i}(1),:),:);
    V = Vandermonde_LP(nDegParametric(i),coord);
    [L,U,P] = lu(V');
    p = orthopoly1D_deriv(pt{i},nDegParametric(i));
    N = U\(L\(P*p));
    for j = 1:nOfTerms
        iterm = reducedBasis{i}(:,j);
        interpolatedParam(j) = N'*iterm(parametricMesh{i,2}(elem{i}(1),:));
    end
    interpolatedCoef = interpolatedCoef.*interpolatedParam;
end

%Interpolate PGD approximation
u = (interpolatedCoef*reducedBasis{end}.').';


function N = nodalConnectivity(T)
 
nNodes = max(max(T));
N = zeros(nNodes,2);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
   Te = T(ielem,:);
   nn_Te = nn(Te);
   for kk = 1:2
       N(Te(kk),nn_Te(kk)) = ielem;
   end
   nn(Te) = nn(Te) + 1;
end


