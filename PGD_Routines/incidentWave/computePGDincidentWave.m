function ui = computePGDincidentWave(PGDmeshes,PGDalgorithm,parameters)

dimTHETA  = findPGDdimension('THETA',parameters.PGDdimensions);
dimK      = findPGDdimension('K',parameters.PGDdimensions);
dimKTHETA = findPGDdimension('KTHETA',parameters.PGDdimensions);

if dimK && dimTHETA
    if parameters.iwparam.svd
        %Use SVD with 2D spatial problems in \Omega_PML x (I_w x I_\theta)
        ui = computeSVDincidentWave_XY_K_THETA(PGDmeshes,PGDalgorithm,parameters);
    else
        %Use PGD projection with 1D spatial problems in {I_x,I_y} x (I_w x I_\theta)
        ui = computePGDincidentWave_XY_K_THETA(PGDmeshes,PGDalgorithm,parameters);
    end
elseif dimKTHETA
    if parameters.iwparam.svd
        %Use SVD with 2D spatial problems in \Omega_PML x (I_w x I_\theta)
        ui = computeSVDincidentWave_XY_KTHETA(PGDmeshes,PGDalgorithm,parameters);
    else
        %Use PGD projection with 1D spatial problems in {I_x,I_y} x (I_w x I_\theta)
        error('PGD Incident wave for parametric dimension KTHETA not implemented yet')
    end
elseif dimK
    if parameters.iwparam.svd
        %Use SVD projection in \Omega_PML x I_w
        ui = computeSVDincidentWave_XY_K(PGDmeshes,parameters);
    else
        %Use PGD projection in \Omega_PML x I_w
        error('PGD Incident wave for parametric dimension K not implemented yet')
    end
elseif dimTHETA
    if parameters.iwparam.svd
        %Use SVD projection in \Omega_PML x I_\theta
        ui = computeSVDincidentWave_XY_THETA(PGDmeshes,parameters);
    else
        %Use PGD projection in \Omega_PML x I_\theta
        error('PGD Incident wave for parametric dimension THETA not implemented yet')
    end
else
    ui = computeIncidentWave_XY(PGDmeshes,parameters);
end



