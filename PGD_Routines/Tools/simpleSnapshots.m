setPGDpath()

% load('midscat_2_KR100_const.mat') %data GUI
load() %PGDmeshes
load() %Parameters

a = 6.2832;
b = 31.4159;
m = 1000;
omegav = linspace(a,b,m);

% ufem = zeros(size(data.mesh.X,1),m);
ufem = zeros(size(PGDmeshes(1).X,1),m);
for i = 1:m
    disp(['-------SNAPSHOT ' num2str(i) ' value ' num2str(omegav(i)) '-------'])
%     data.bottom = rmfield(data.bottom,'ccg');
%     data.ip.period = 2*pi/omegav(i);
%     data = run_domain(data,[]);
%     data = run_boundary(data,[]);
%     ufem(:,i) = data.solution;
    computeFEMfromPGDdata([omegav(i),parameters.fixedParametricDims(2)],PGDmeshes,parameters);
end

bigsave('ufemdata','ufem')