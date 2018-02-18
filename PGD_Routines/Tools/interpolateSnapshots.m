clear all
setPGDpath()

%% Loads

load('mataro[L200][T10to16][N8][P4]_const.mat')
ufem_e = load('ufemdataMataro_10x10_ufem.mat'); %at interpolation points
% ufem = load('ufemdataMataro_25x25.mat'); %at sample points
load('PGDdata.mataro.K_THETA[10to16][190to270][50x50]_PGDmeshes.mat') %PGD meshes
load('RB.mataro.K_THETA[10to16][190to270][50x50]_i3_pg_pgd.mat') %PGD reduced basis

%% User data

%Solver
solver = 'pgd';
method = 'linear'; %(if solver == 'interp')

nOfPGDtermsVEC = size(pgd.RB{1},2); %(if solver == 'pgd', empty otherwise)
if nOfPGDtermsVEC(end) < size(pgd.RB{1},2), nOfPGDtermsVEC = [nOfPGDtermsVEC,size(pgd.RB{1},2)]; end

DoPlots = true;
outputFile = 'RB.mataro.K_THETA[10to16][190to270][50x50]_i3_pg[p1].mat';

%Limit error for plots
lim = 0.05;

%Error to be computed in the full parametric domain
errorType = 1; %1 for H(x = a) ; 2 for H_max(x\in T) ; 3 for norm H(x\in T)
dataset = {
%             5483    %punto del espigon de entrada de mataro
            5509    %punto delante de la bocana de mataro
%             12585   %punto interior profundo de mataro
%             10527   %punto interior primero de mataro
%             'T_1'   %zona delante de la bocana de mataro
%             'T_2'   %zona interior primera de mataro
%             'T_3'   %zona interior profunda de mataro
           };
       
%Parametric domain: limits for frequency (1) and direction (2)
a1 = 0.3927;
b1 = 0.6283;
a2 = 3.3161;
b2 = 4.7124;

%% Interpolation procedure

%dataset info
if errorType == 1
    Tnodes = dataset{1};
elseif any(errorType == [2,3])
    Tname = dataset{1};
    Tnodes = unique(data.mesh.(Tname));
end
nTnodes = numel(Tnodes);

%Exact FEM solution at interpolation points
m_e1 = size(ufem_e.ufem,2);
m_e2 = size(ufem_e.ufem,3);
v1i = linspace(a1,b1,m_e1);
v2i = linspace(a2,b2,m_e2);
[Xi,Yi] = meshgrid(v1i,v2i);
sol_e = abs(ufem_e.ufem);
if errorType == 1
    sol_e = reshape(sol_e(Tnodes,:,:),m_e1,m_e2);
elseif errorType == 2
    sol_e = max(sol_e(Tnodes,:,:),[],1);
    sol_e = reshape(sol_e,m_e1,m_e2);
elseif errorType == 3
    sol_e = reshape(sol_e(Tnodes,:,:),nTnodes,m_e1,m_e2);
end

%Interpolates solutions (solver-dependent)
for j = 1:numel(nOfPGDtermsVEC)
    nOfPGDterms = nOfPGDtermsVEC(j);
    disp(['Number of terms = ' num2str(nOfPGDterms)])
    
    if strcmpi(solver,'interp') %independent from nOfPGDterms

        %FEM solution at sample points
        m1 = size(ufem.ufem,2);
        m2 = size(ufem.ufem,3);
        v1 = linspace(a1,b1,m1);
        v2 = linspace(a2,b2,m2);
        [X,Y] = meshgrid(v1,v2);
        sol = abs(ufem.ufem); 
        if errorType == 1
            sol = reshape(sol(Tnodes,:,:),m1,m2);
        elseif any(errorType == [2,3])
            sol = reshape(sol(Tnodes,:,:),nTnodes,m1,m2);
        end

        %Interpolated FEM solution
        if errorType == 1
            Zi = interp2(X,Y,sol,Xi,Yi,method);
        elseif any(errorType == [2,3])
            Zi = zeros(nTnodes,m_e1,m_e2);
            for i = 1:nTnodes
                Zi(i,:,:) = interp2(X,Y,reshape(sol(i,:,:),m1,m2),Xi,Yi,method);
            end
            if errorType == 2
                Zi = reshape(max(Zi,[],1),m_e1,m_e2);
            end
        end
        
    elseif strcmpi(solver,'pgd')
        if exist('projectedpgd','var'), pgd = projectedpgd; clear('projectedpgd'), end

        %Interpolated PGD solution  u = interpolatePGD(dim=3.5,param1,param2,S,R,X_s,T_s,Deg_s,Xp,Np)
        if exist('PROJmeshes','var')
            error('not implemented yet')
%                 (3.5,v1i',v2i',...
%                 pgd.RB{2},pgd.RB{1},...
%                 PROJmeshes.X,PROJmeshes.T,PROJmeshes.referenceElement.degree,...
%                 Tnodes,nOfPGDterms);
        else
            Xit = Xi';
            Yit = Yi';
            Zi = interpolatePGD(pgd,[2 3],{Xit(:) Yit(:)},PGDmeshes,nOfPGDterms);
        end
        Zi = Zi(Tnodes,:);
        Zi = abs(Zi);

        if errorType == 1
            Zi = reshape(Zi,m_e1,m_e2);
        elseif errorType == 2
            Zi = reshape(max(Zi,[],1),m_e1,m_e2);
        elseif errorType == 3
            Zi = reshape(Zi,nTnodes,m_e1,m_e2);
        end
    end

    %% Errors

    %Error
    e = abs(Zi - sol_e) ./ abs(sol_e);
    e_norm = norm(Zi(:) - sol_e(:)) / norm(sol_e(:));
    e_norm %#ok<*NOPTS>

    %Binary error
    ebounded = e;
    pos1 = e > lim; ebounded(pos1) = 0;
    ebounded(~pos1) = NaN;
    perc = 100 * numel(find(pos1 == true)) / numel(e);
    perc

    %% Plots

    if DoPlots && errorType ~= 3 
        figure, h1 = surf(v1i,v2i,e); view(2)
        axis tight
        caxis_lims = [0 lim];
        caxis(caxis_lims)

        figure, h2 = surf(v1i,v2i,ebounded); view(2), grid('off'), box('on')
        colormap([1 1 1; 0 0 0])
        axis tight
    end
    
    %Save error
    errorVEC(j) = e_norm;
    nVEC = nOfPGDtermsVEC(1:j);
    if ~isempty(outputFile), save(outputFile,'errorVEC','nVEC'), end
end





