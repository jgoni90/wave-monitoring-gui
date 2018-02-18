%SCRIPT TO SETUP THE PGD CODE: dimensions, meshes, parameters, algorithm...

%% SETUP NEW PGD OR LOAD / SAVE PGDdata FROM FILE

%Load flag oder: [coefficients, matrices and RHS]
loadPGDfromFile = [0,0];

%Save flag oder: [meshes, coefficients, matrices and RHS]
%matrices and RHS == 2 only saves the mass matrices and ingnores the RHS
savePGDintoFile = [1,1,2];

%File name
PGDfile = 'PGDdata.mataro.K_THETA[30to333.3][190to270][100x100]';

%% PGD PARAMETERS

%General parameters
%ResidualType LMC  = Last Mode Contribution 
%ResidualType DRFS = Discretized Residual in the Full Space
maxterms            = 2500;
maxiter             = 3;
toliter             = 1e-3;
tolu                = 1e-8;
residualType        = 'LMC';
residualEachTerm    = 1;
saveEachTerm        = 25;
outputFile          = 'RB.mataro.K_THETA[30to333.3][190to270][100x100]_pg_i3';
recoverPGD          = false;
recoverFile         = '';

%Parameters for the variable coefficients
coefparam.svd        = true;
coefparam.load       = 'PGD_freq_depth_T3_T16_d1_d65_1000x1000.mat';
coefparam.relerror   = 1e-5;

%Parameters for the incident wave
iwparam.load         = '';
iwparam.save         = '';
iwparam.svd          = true;
iwparam.loadfileX    = 'PGD_UIX_X_KTHETA_1000_500x500[0to1][1to600][PIto2PI].mat';
iwparam.loadfileY    = 'PGD_UIY_Y_KTHETA_1000_500x500[0to1][1to600][PIto2PI].mat';
iwparam.maxterms     = 100;
iwparam.maxiter      = 5;
iwparam.resisualType = 'ALPHACOEF';
iwparam.relerror     = 1e-5;
iwparam.maxGBmemory  = 7;

%% PGD TYPE

%ALPHA_(boundaryName) dimensions have the same boundary order as in fileloadmeshes{1}

PGDdimensions = {    %Name    %Dim  %Load with file  %parametric dim
%                      'X'        1        0                0          
                     'XY'       2        1                0       
                     'K'        1        0                1
                     'THETA'    1        0                1
%                      'KTHETA'   2        0                1
%                      'ALPHA_one'        1        0                1
%                      'ALPHA_two'        1        0                1
%                      'ALPHA_three'        1        0                1
%                      'ALPHA_four'        1        0                1
                 };

%% MESH FILES

%Non-empty is mandatory (Berkhoff GUI file) for XY dimension
%Empty is optional for the rest of dimensions

fileloadmeshes = {
                  'mataro[L200][T10to16][N8][P4]_const.mat'
                  []
                  []
%                   []
%                   []
                 };
             
%% PGD ALGORITHM

%Petrov-Galerkin
PGDalgorithm.PG                         = true;

%Dual
PGDalgorithm.dual.value                 = false;
PGDalgorithm.dual.RBfile                = 'RB_mataro4D_i25_pg.mat';
PGDalgorithm.dual.type                  = 'area'; %line or area
PGDalgorithm.dual.spatialAttributeName  = 'T_errorDual'; 

%Update
PGDalgorithm.update.value               = false;
PGDalgorithm.update.parameter           = {'K'};
PGDalgorithm.update.parameterStoreRB    = false;

%L2-Projection
PGDalgorithm.projection.value           = false;
PGDalgorithm.projection.projectEachTerm = 10;
PGDalgorithm.projection.maxiter         = 30;
PGDalgorithm.projection.residualType    = 'LMC';
PGDalgorithm.projection.projCoeff       = false;
PGDalgorithm.projection.relerror        = 1e-6;
PGDalgorithm.projection.tolfactor       = 1.02;

%% MESH PARAMETERS

%For loaded meshes:     empty
%For 1D X mesh:         Xini     Xend      nOfnodes nDeg nDegNodes
%For 1D K mesh:         Tini     Tend      nOfNodes nDeg nDegNodes
%For 1D THETA mesh:     THETAini THETAend  nOfNodes nDeg nDegNodes
%For 2D KTHETA mesh:    Tini Tend nk THETAini THETAend nTHETA nDeg nDegNodes equalSpacedNodes
%For 1D ALPHA mesh:     ALPHAini ALPHAend  nOfNodes nDeg nDegNodes

meshparameters = {
                  {[]}
                  {30    333.3     100    1    3}
                  {190      270   100    1    3}
%                   {0.4    1    19    0    360    19    1    4    false}
%                   {0    1     10    1    3}
%                   {0    1     10    1    3}
%                   {0    1     10    1    3}
%                   {0    1     10    1    3}
                 };
 
%Vector of fixed parametric dims, in order: [K THETA ALPHA1 ... ALPHAn]
%K in [1/s], THETA in [rad] and ALPHA1,...,ALPHAn with the same order as in fileloadmeshes{1}
%-1 value for those dims which are variable
%ALPHA1 = -2 indicates that all the ALPHA fixed values are obtained from fileloadmeshes{1} 

fixedParametricDims = [-1 -1 -2];

%% SETUP PGD

run setupPGD

%% SAVE

if any(savePGDintoFile)
    disp('SAVING DATA STRUCTURE')
    savestring        = {'PGDmeshes','PGDcoefs','PGDmatrices'    ,'PGDrhs'};
    savestringStorage = {''         ,''        ,'PGDmassmatrices',''};
    savePGDintoFile = [savePGDintoFile(1:3),savePGDintoFile(3),savePGDintoFile(4:end)];
    bigsave(PGDfile,savestring{savePGDintoFile == 1})
    bigsave(PGDfile,savestringStorage{savePGDintoFile == 2})
end
disp('DONE!')
    
    
    
