function ui = computePGDcoefs(parameters,PGDmeshes)

dimK = findPGDdimension('K',parameters.PGDdimensions);

if dimK
    if parameters.coefparam.svd
        %Use SVD in \Omega x I_w
        ui = computeSVDcoefs_XY_K(parameters,PGDmeshes);
    else
        %Use PGD projection
        error('PGD coefs for parametric dimension K not implemented yet')
    end
else 
    ui = computeCoefs_XY(parameters);
end