function fi = projectQUADonQUADmesh(f,SAMPLEPOINTS,X2,T2,degree)

% fi = projectQUADonQUADmesh(f,sp,X,T,re)
%
% Interpolates a function f(x,y) given on a rectangular 2D QUAD mesh defined 
% by X, T and re (referenceElement) over a 2D rectangular mesh given by 1D 
% vectors x = sp{1} and y = sp{2}; f is the vector or matrix of nodal values. 
% The mesh should be sorted in row-ascend order from (0,0) to (end,end).
% Result will be returned with the same row-ascend order.

%Initialize
nsd             = 2; %spatial dimensions

osp             = cell(nsd,1);
nOfsp           = zeros(nsd,1);
for j = 1:nsd
    osp{j}      = sort(SAMPLEPOINTS{j},'ascend');
    nOfsp(j)    = size(osp{j},1);
end
poselem         = cell(nsd,1);
elems           = cell(nsd,1);
xisp            = cell(nsd,1);
nnodesQUAD      = size(T2,2);
nOffun          = size(f,2);
fi              = zeros(prod(nOfsp),nOffun);
oneDimMap       = @(d1,d2,x)(1/(d1(1)-d1(2)))*((d2(1)-d2(2))*x + d2(2)*d1(1) - d2(1)*d1(2));

%1D vectors get from the QUAD mesh (row-ascend order)
X               = cell(nsd,1);
T               = cell(nsd,1);
nnodes          = zeros(nsd,1);
nelems          = zeros(nsd,1);
for i = 1:nsd
    X{i}        = unique(X2(:,i));
    nnodes(i)   = size(X{i},1);
    nelems(i)   = (nnodes(i)-1) / degree;
    T{i}        = create1Dconec(nelems(i),degree);
end
nen1D           = size(T{1},2);
N               = cell(nsd,1);
NQUAD           = zeros(nOfsp(1),nOfsp(2),nnodesQUAD);

switch degree
    case 1
        
        %Search for 1D elements and compute samplepoints on parametric space
        for j = 1:nsd
            
            %Initialization
            poselem{j}    = zeros(nelems(j),2);
            elems{j}      = zeros(nOfsp(j),1);
            xisp{j}       = zeros(nOfsp(j),1);
            pos           = 1;
            prepos        = 1;
            tolv          = zeros(nelems(j),1);
            tolv(end)     = 1e-8;
            
            %Loop in 1D elements
            for i = 1:nelems(j)
                Xe = X{j}(T{j}(i,:));
                x2 = Xe(2);
                tol = tolv(i);
                while pos <= nOfsp(j) && osp{j}(pos) < x2+tol, pos = pos + 1; end
                poselem{j}(i,1) = prepos;
                poselem{j}(i,2) = pos-1;
                posv = prepos:pos-1;
                elems{j}(posv) = i;
                xisp{j}(posv) = oneDimMap(Xe,[-1,1],osp{j}(posv));
                prepos = pos;
            end
            
            %Computes 1D linear shape functions on [-1,1] at samplePoints
            N{j} = (1/2)*[1-xisp{j}, 1+xisp{j}];
        end
        
        %QUAD shape functions in the reference element
        tensorN = bsxfun(@times,reshape(N{1},nOfsp(1),1,nen1D,1),...
            reshape(N{2},1,nOfsp(2),1,nen1D));
        NQUAD(:,:,1) = tensorN(:,:,1,1);
        NQUAD(:,:,2) = tensorN(:,:,2,1);
        NQUAD(:,:,3) = tensorN(:,:,2,2);
        NQUAD(:,:,4) = tensorN(:,:,1,2);
        
        %Loop in X elements
        uelems1 = unique(elems{1})';
        uelems2 = unique(elems{2})';
        for i = uelems1
            xini = poselem{1}(i,1);
            xfin = poselem{1}(i,2);
            xsample = (xini:xfin)';
            nOfxsample = xfin-xini+1;
            onesx = ones(nOfxsample,1);
            
            %Loop in Y elements
            for j = uelems2
                yini = poselem{2}(j,1);
                yfin = poselem{2}(j,2);
                ysample = yini:yfin;
                nOfysample = yfin-yini+1;
                onesy = ones(1,nOfysample);
                
                %Interpolation
                QUADelem = i + (j-1)*nelems(1);
                jnodes = T2(QUADelem,:);
                jf = reshape(NQUAD(xsample,ysample,:),nOfxsample*nOfysample,nnodesQUAD) * f(jnodes,:);
             
                %Storage
                A_xsample = xsample(:,onesy);
                A_xsample = A_xsample(:);
                A_ysample = ysample(onesx,:);
                A_ysample = A_ysample(:);
                samplepos = A_xsample + (A_ysample-1)*nOfsp(1);
                fi(samplepos,:) = jf;
            end
        end
        
    otherwise
        
        error('Interpolation for degree > 1 not implemented yet')
end





