function u = interpolatePGD(dim,varargin)

% Return u(x,y) = u(x,y,param1,param2), interpolated at parameter values param1 
% and param2, from reduced basis R(x,y)*S(param1)*T(param2).
%
%       u = interpolatePGD(dim=4,param1,param2,S,T,R,X_s,T_s,X_t,T_t,Xp,Np)
%       u = interpolatePGD(dim=3.5,param1,param2,S,R,X_s,T_s,Deg_s,Xp,Np)
%       u = interpolatePGD(dim=3,param1,S,R,X_s,T_s,Xp,Np)
%
%   Inputs:
%       dim: number of total dimensions (space + parameters).
%       param: parameter vector where the PGD approximation will be
%              interpolated.
%       S,T,R: reduced basis (R(x,y) has to be the last one).
%       X,T: one-dimensional parametric mesh, X: nodal position, T: connectivity.
%       Xp: list of spatial nodes (empty for Xp = X).
%       Np: number of terms of the reduced basis (empty for Np = size(reduced basis,2))

%Check dimension type
switch dim
    case 3
        pt{1} = varargin{1};
        reducedBasis(1) = varargin(2);
        reducedBasis(2) = varargin(3);
        parametricMesh(1,1) = varargin(4);
        parametricMesh(1,2) = varargin(5);
        nOfParametricDimensions = 1;
    case 3.5
        pt{1} = varargin{1};
        pt{2} = varargin{2};
        reducedBasis(1) = varargin(3);
        reducedBasis(2) = varargin(4);
        parametricMesh(1,1) = varargin(5);
        parametricMesh(1,2) = varargin(6);
        deg = varargin{7};
        nOfParametricDimensions = 1;
    case 4
        pt{1} = varargin{1};
        pt{2} = varargin{2};
        reducedBasis(1) = varargin(3);
        reducedBasis(2) = varargin(4);
        reducedBasis(3) = varargin(5);
        parametricMesh(1,1) = varargin(6);
        parametricMesh(1,2) = varargin(7);
        parametricMesh(2,1) = varargin(8);
        parametricMesh(2,2) = varargin(9);
        nOfParametricDimensions = 2;
end
xnodes = varargin{end-1};
nterms = varargin{end};
if isempty(nterms), nterms = size(reducedBasis{1},2); end

%Interpolate non-spacial dimensions
interpolatedParam = cell(nOfParametricDimensions,1);
if dim == 3.5
    for i = 1:nOfParametricDimensions
        interpolatedParam{i} = projectQUADonQUADmesh(reducedBasis{i}(:,1:nterms),...
            {pt{1},pt{2}},parametricMesh{i,1},parametricMesh{i,2},deg);
    end
else
    for i = 1:nOfParametricDimensions
        deg = size(parametricMesh{i,2},2) - 1;
        interpolatedParam{i} = interpolateOn1Dmesh(reducedBasis{i}(:,1:nterms),pt{i},...
            parametricMesh{i,1},parametricMesh{i,2},...
            deg);
    end
end

%Interpolated reduced basis u(x,y) for parametric points pt{1} x ... x pt{n}
if isempty(xnodes), RB = reducedBasis{end}(:,1:nterms).'; else RB = reducedBasis{end}(xnodes,1:nterms).'; end
nx = size(RB,2);
if nOfParametricDimensions == 1
    np1 = size(interpolatedParam{1},1);
    u = zeros(nx,np1);
    for j = 1:np1
        u(:,j) = (interpolatedParam{1}(j,:)*RB).';
    end
elseif nOfParametricDimensions == 2
    np1 = size(interpolatedParam{1},1);
    np2 = size(interpolatedParam{2},1);
    u = zeros(nx,np1,np2);
    for j = 1:np1
        for k = 1:np2
            coef = interpolatedParam{1}(j,:).*interpolatedParam{2}(k,:);
            u(:,j,k) = (coef*RB).';
        end
    end
end


