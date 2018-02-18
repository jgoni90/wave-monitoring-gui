function ui = computeIncidentWave_XY(PGDmeshes,parameters)

%% INFORMATION OF THE INCIDENT WAVE FUNCTION

omega = parameters.fixedParametricDims(1);
theta = parameters.fixedParametricDims(2);
dimXY = findPGDdimension('XY',parameters.PGDdimensions);

info.PGDmeshesField = {'T','ext'};

%% COMPUTATION OF THE INCIDENT WAVE

T = myGetField(PGDmeshes(dimXY),info.PGDmeshesField);
xy_nodes = unique(T);
ui.meshes(dimXY).nodes = xy_nodes;
bottom = parameters.meshes(dimXY).bottom(xy_nodes(1));
k = computeWaveNumber(omega,bottom);

ui.nOfTerms = 1;
x = PGDmeshes(dimXY).X(xy_nodes,1);
y = PGDmeshes(dimXY).X(xy_nodes,2);
ui.RB{dimXY} = exp(1i*k*(x*cos(theta) + y*sin(theta)));

%% OTHER FUNCTIONS

function v = myGetField(s,f), v = s; for i = 1:numel(f), v = v.(f{i}); end



