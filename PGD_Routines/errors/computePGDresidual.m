function [res,resTensor] = computePGDresidual(parameters,pgd,meshes,matrices,rhs,algorithm,flag)

%% COMPUTE THE WEAK RESIDUAL

Uterm = pgd.counters.Uterm;
parametricDims = findPGDdimension('PARAMETRICDIM',parameters.PGDdimensions);
nOfParametricDims = length(parametricDims);

dimX = findPGDdimension('X',parameters.PGDdimensions);
dimXY = findPGDdimension('XY',parameters.PGDdimensions);

im = sqrt(-1);

%----

if dimX
    
end

%----

if dimXY
    
    if flag
        iniTerm = pgd.counters.residualUterm;
        prevresMat = pgd.errors.residualMat;
    else
        iniTerm = 1;
        prevresMat = 0;
    end
    
    switch nOfParametricDims
        
        %-----
        
        case 1
            
            dim = parameters.PGDdimensions{parametricDims,1};
            sizeTensor = size(meshes(dimXY).X,1)*size(meshes(parametricDims).X,1);
            
            if any(strcmp(dim,{'K' 'KTHETA'}))
                
                if algorithm.dual.value
                    resMat = 0;
                    for i = iniTerm:Uterm
                        resMat = resMat + matrices(dimXY).K*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).M*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat - matrices(dimXY).M*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mk2*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat + im*matrices(dimXY).C*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mk*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat - matrices(dimXY).B*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mkr*pgd.RB{parametricDims}(:,i)).';
                    end
                    resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                    clear resMat prevresMat
                    rhsTensor = reshape(rhs,sizeTensor,1);
                    res = norm(resTensor - rhsTensor) / norm(rhsTensor);
                    
                else
                    resMat = 0;
                    for i = iniTerm:Uterm
                        resMat = resMat + matrices(dimXY).K*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).M*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat - matrices(dimXY).M*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mk2*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat - im*matrices(dimXY).C*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mk*pgd.RB{parametricDims}(:,i)).';
                        resMat = resMat - matrices(dimXY).B*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).Mkr*pgd.RB{parametricDims}(:,i)).';
                    end
                    resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                    clear resMat prevresMat
                    rhsTensor = reshape(rhs,sizeTensor,1);
                    res = norm(resTensor + rhsTensor) / norm(rhsTensor);
                end
                
            elseif strcmp(dim,'THETA')
                
                k = parameters.fixedParametricDims(1);
                if algorithm.dual.value
                    paramrad = (-im*k - 1/(2*parameters.rad));
                    A_xy = (matrices(dimXY).K - k^2*matrices(dimXY).M + ...
                            im*k*matrices(dimXY).C - paramrad*matrices(dimXY).B);
                        
                    resMat = 0;
                    for i = iniTerm:Uterm
                        resMat = resMat + A_xy*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).M*pgd.RB{parametricDims}(:,i)).';
                    end
                    resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                    clear resMat prevresMat
                    rhsTensor = reshape(rhs,sizeTensor,1);
                    res = norm(resTensor - rhsTensor) / norm(rhsTensor);
                    
                else
                    paramrad = (im*k - 1/(2*parameters.rad));
                    A_xy = (matrices(dimXY).K - k^2*matrices(dimXY).M - ...
                            im*k*matrices(dimXY).C - paramrad*matrices(dimXY).B);

                    resMat = 0;
                    for i = iniTerm:Uterm
                        resMat = resMat + A_xy*pgd.RB{dimXY}(:,i) *...
                                  (matrices(parametricDims).M*pgd.RB{parametricDims}(:,i)).';
                    end
                    resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                    clear resMat prevresMat
                    rhsTensor = reshape(rhs,sizeTensor,1);
                    res = norm(resTensor + rhsTensor) / norm(rhsTensor);
                end
                
            end
            
        %-----
            
        case 2
            
            dimK        = findPGDdimension('K',parameters.PGDdimensions);
            dimTHETA    = findPGDdimension('THETA',parameters.PGDdimensions);
            nNodesXY    = size(meshes(dimXY).X,1);
            nNodesK     = size(meshes(dimK).X,1);
            nNodesTHETA = size(meshes(dimTHETA).X,1);
            sizeTensor  = nNodesXY*nNodesK*nNodesTHETA;
            
            if algorithm.dual.value
                resMat = 0;
                for i = iniTerm:Uterm
                    v_xy    = computePGDrightVecs(pgd.RB{dimXY}(:,i),matrices(dimXY));
                    v_k     = computePGDrightVecs(pgd.RB{dimK}(:,i),matrices(dimK));
                    v_theta = computePGDrightVecs(pgd.RB{dimTHETA}(:,i),matrices(dimTHETA));

                    resMat0 = v_xy.K*(v_k.M).' - v_xy.M*(v_k.Mk2).' + im*v_xy.C*(v_k.Mk).' -...
                        v_xy.B*(v_k.Mkr).';
                    resMat = resMat + (v_theta.M)*reshape(resMat0,1,nNodesXY*nNodesK);
                end
                resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                clear resMat prevresMat
                rhsTensor = reshape(rhs,sizeTensor,1);
                res = norm(resTensor - rhsTensor) / norm(rhsTensor);
                
            else
                resMat = 0;
                for i = iniTerm:Uterm
                    v_xy    = computePGDrightVecs(pgd.RB{dimXY}(:,i),matrices(dimXY));
                    v_k     = computePGDrightVecs(pgd.RB{dimK}(:,i),matrices(dimK));
                    v_theta = computePGDrightVecs(pgd.RB{dimTHETA}(:,i),matrices(dimTHETA));

                    resMat0 = v_xy.K*(v_k.M).' - v_xy.M*(v_k.Mk2).' - im*v_xy.C*(v_k.Mk).' -...
                        v_xy.B*(v_k.Mkr).';
                    resMat = resMat + (v_theta.M)*reshape(resMat0,1,nNodesXY*nNodesK);
                end
                resTensor = reshape(resMat,sizeTensor,1) + prevresMat;
                clear resMat prevresMat
                rhsTensor = reshape(permute(rhs,[1 3 2]),sizeTensor,1);
                res = norm(resTensor + rhsTensor) / norm(rhsTensor);
            end
            
    end
    
end
