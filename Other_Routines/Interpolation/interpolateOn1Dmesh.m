function [fi,structure] = interpolateOn1Dmesh(f,SAMPLEPOINTS,X,T,degree,structure)

% fi = interpolateOn1Dmesh(f,samplePoints,X,T,re)
%
% Interpolates a 1D function f(x) at x = samplePoints over a 1D mesh
% defined by X, T and re (referenceElement); f is the vector of nodal values.
% The mesh should be sorted in ascend order.

nOffun      = size(f,2);
nOfsp       = numel(SAMPLEPOINTS);
[osp,opos]  = sort(SAMPLEPOINTS,'ascend');
nElems      = size(T,1);
elems       = zeros(1,nOfsp);
xisp        = zeros(nOfsp,1);
poselem     = zeros(nElems,2);
fi          = zeros(nOfsp,nOffun);
oneDimMap   = @(d1,d2,x)(1/(d1(1)-d1(2)))*((d2(1)-d2(2))*x + d2(2)*d1(1) - d2(1)*d1(2));

switch degree
    case 1
        
        %Search for 1D elements and compute samplepoints on parametric space
        if nargin == 5 %do computation
            pos = 1;
            prepos = 1;
            tolv = zeros(nElems,1);
            tolv(end) = 1e-4;
            for i = 1:nElems
                Xe = X(T(i,:));
                x2 = Xe(2);
                tol = tolv(i);
                while pos <= nOfsp && osp(pos) <= x2+tol, pos = pos + 1; end
                poselem(i,1) = prepos;
                poselem(i,2) = pos-1;
                posv = prepos:pos-1;
                elems(posv) = i;
                xisp(posv) = oneDimMap(Xe,[-1,1],osp(posv));
                prepos = pos;
            end
            
            if nargout == 2 %use the function to export basic data for further interpolation
                structure.elems = elems;
                structure.poselem = poselem;
                return
            end    
        elseif nargin == 6 %load basic data for interpolation
            elems = structure.elems;
            poselem = structure.poselem;
        end
        
        %Computes 1D linear shape functions on [-1,1] at samplePoints
        N = (1/2)*[1-xisp, 1+xisp];
        
        %Interpolate
        for j = unique(elems)
            jnodes = T(j,:);
            jsample = poselem(j,1):poselem(j,2);
            fi(opos(jsample),:) = N(jsample,:) * f(jnodes,:);
        end
        
    otherwise
        
        error('Interpolation for degree > 1 not implemented yet')
end


