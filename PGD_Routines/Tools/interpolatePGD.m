function u = interpolatePGD(pgd,dims,pts,meshes,nterms)

% u = interpolatePGD(pgd,dims,pts,meshes,nterms) returs u(:,i) = pgd(x,y;pts{i}(1),...,pts{i}(n)).
% INPUT:
%     pgd:    pgd structure
%     dims:   parametric dimensions (non-spatial 1D or 2D) to be interpolated
%     pts:    parametric points where the pgd should be interpolated
%     meshes: PGDmeshes structure
%     nterms: (OPTIONAL) number of PGD terms to be used. DEFAULT: all
%
% Limitations:
%     - Linear/bilinear parametric meshes.
%
% Example: Consider u(x,y,p1,p2,p3) = sum_{i=1}^N F1(p1,p2)*F2(x,y)*F3(p3)
%          
%          Case 1: interpolate at p1 = 2, p2 = 3, p3 = 6
%                  TYPE: u = interpolatePGD(pgd,[1,3],{[2,3],6},meshes)
%          Case 2: interpolate at p1 = [2,1], p2 = [3,7], p3 = [6,9]
%                  TYPE: u = interpolatePGD(pgd,[1,3],{[2,1;3,7],[6,9]},meshes)

%Sizes
ndims   = length(dims);
npoints = length(pts{1}(:,1));
if nargin < 5, nterms = size(pgd.RB{1},2); end

%Interpolate parametric dimensions
coef = ones(npoints,nterms);
for i = 1:ndims
    idim = dims(i);
    deg = meshes(idim).referenceElement.degree;
    
    %1D or 2D interpolation
    if size(meshes(idim).X,2) == 1 %1D parametric mesh
        icoef = interpolateOn1Dmesh(pgd.RB{idim}(:,1:nterms),pts{i},...
            meshes(idim).X,meshes(idim).T,deg);
    elseif size(meshes(idim).X,2) == 2 %2D parametric mesh
        icoef = projectQUADonQUADmesh(pgd.RB{idim}(:,1:nterms),...
            {pts{i}(:,1),pts{i}(:,2)},meshes(idim).X,meshes(idim).T,deg);
        icoef = icoef(1:npoints+1:end,:); %get diagonal values
    else 
        error('Only 1D and 2D parametric meshes are allowed')
    end
    
    %Update coefficient
    coef = coef .* icoef;
end

%Interpolated PGD
dimXY = setdiff(1:length(pgd.RB),dims);
u = pgd.RB{dimXY}(:,1:nterms)*coef.';