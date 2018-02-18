
%% OPTIONS

whichVariableDoIwant2interpolate    = 'K';
paramValues                         = 6.32;
runFEM                              = 1;
intervalPGDterms                    = 1;
maxTerms                            = 9;
PGDfile                             = '';
PGDdatafile                         = '';

%% DATA

if ~isempty(PGDfile),       load(PGDfile),      end
if ~isempty(PGDdatafile),   load(PGDdatafile),  end
if exist('PGDmassmatrices','var'), PGDmatrices = PGDmassmatrices; end
if ~isempty(maxTerms), for i = 1:parameters.nOfPGDdimensions, pgd.RB{i} = pgd.RB{i}(:,1:maxTerms); end, end

%% CODE

dimXY = findPGDdimension('XY',parameters.PGDdimensions);
switch whichVariableDoIwant2interpolate
    
    case 'K' %%%% INTERPOLATION FOR FREQUECY
        dimK = findPGDdimension('K',parameters.PGDdimensions);

            %FEM
            if runFEM
                snapshotsValue = [paramValues,parameters.fixedParametricDims(2)];
                usol = computeFEMfromPGDdata(snapshotsValue,PGDmeshes,parameters);
            end
            uzeros = zeros(size(usol));
            nodesInt = unique(PGDmeshes(dimXY).T.int);
            ufem = uzeros;
            %ufem(nodesInt) = usol(nodesInt);
            ufem(nodesInt) = abs(usol(nodesInt));
            normcoef = sqrt(ufem' * PGDmatrices(dimXY).M * ufem);

            %PGD
            PGDterms  = size(pgd.RB{1},2);
            PGDtermsv = 1:intervalPGDterms:PGDterms;
            relerr    = zeros(length(PGDtermsv),1);
            for i = 1:length(PGDtermsv)
                disp(['Interpolating with ' num2str(PGDtermsv(i)) ' terms...'])
                pos = 1:PGDtermsv(i);
                uint = uzeros;
                u = interpolatePGD...
                    (3,paramValues,...
                    pgd.RB{dimK}(:,pos),pgd.RB{dimXY}(:,pos),...
                    PGDmeshes(dimK).X,PGDmeshes(dimK).T,PGDmeshes(dimK).referenceElement.degree,...
                    [],PGDtermsv(i));
                uint(nodesInt) = u(nodesInt);
                %udiff = uint-ufem;
                udiff = abs(uint) - ufem;
                relerr(i) = sqrt(udiff' * PGDmatrices(dimXY).M * udiff) / normcoef;
                %relerr(i) = abs((max(abs(uint))-max(ufem)) / max(ufem));
            end
            
     case 'K_THETA' %%%% INTERPOLATION FOR FREQUECY AND ANGLE OF INCIDENCE
        dimK = findPGDdimension('K',parameters.PGDdimensions);
        dimTHETA = findPGDdimension('THETA',parameters.PGDdimensions);
        
            %FEM
            if runFEM
                snapshotsValue = [paramValues(1),paramValues(2)];
                usol = computeFEMfromPGDdata(snapshotsValue,PGDmeshes,parameters);
            end
            uzeros = zeros(size(usol));
            nodesInt = unique(PGDmeshes(dimXY).T.int);
            ufem = uzeros;
            %ufem(nodesInt) = usol(nodesInt);
            ufem(nodesInt) = abs(usol(nodesInt));
            normcoef = sqrt(ufem' * PGDmatrices(dimXY).M * ufem);
            
            %PGD
            PGDterms  = size(pgd.RB{1},2);
            PGDtermsv = 1:intervalPGDterms:PGDterms;
            relerr    = zeros(length(PGDtermsv),1);
            for i = 1:length(PGDtermsv)
                disp(['Interpolating with ' num2str(PGDtermsv(i)) ' terms...'])
                pos = 1:PGDtermsv(i);
                uint = uzeros;
                u = interpolatePGD...
                    (4,paramValues(1),paramValues(2),...
                    pgd.RB{dimK}(:,pos),pgd.RB{dimTHETA}(:,pos),pgd.RB{dimXY}(:,pos),...
                    PGDmeshes(dimK).X,PGDmeshes(dimK).T,PGDmeshes(dimK).referenceElement.degree,...
                    PGDmeshes(dimTHETA).X,PGDmeshes(dimTHETA).T,PGDmeshes(dimTHETA).referenceElement.degree,...
                    [],PGDtermsv(i));
                uint(nodesInt) = u(nodesInt);
                %udiff = uint-ufem;
                udiff = abs(uint) - ufem;
                relerr(i) = sqrt(udiff' * PGDmatrices(dimXY).M * udiff) / normcoef;
                %relerr(i) = abs((max(abs(uint))-max(ufem)) / max(ufem));
            end
end

%Plot PGD solution
h = semilogy(PGDtermsv,relerr);
hold on

