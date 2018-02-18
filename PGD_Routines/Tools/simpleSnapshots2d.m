setPGDpath()

% load('mataro[L200][T10to16][N8][P4]_const.mat') %Data GUI
load('PGDdata.mataro.K_THETA[10to16][190to270][50x50]_PGDmeshes.mat') %PGD meshes
load('RB.mataro.K_THETA[10to16][190to270][50x50]_i3_pg_parameters.mat') %Parameters

a1 = 0.3927;
b1 = 0.6283;
m1 = 10;
a2 = 3.3161;
b2 = 4.7124;
m2 = 10;
v1 = linspace(a1,b1,m1);
v2 = linspace(a2,b2,m2);

% ufem = zeros(size(data.mesh.X,1),m1,m2);
ufem = zeros(size(PGDmeshes(1).X,1),m1,m2);
cont = 1;
for i = 1:m1
%     data.ip.period = 2*pi/v1(i);
%     data.bottom = rmfield(data.bottom,'ccg');
    for j = 1:m2
        disp(['-------SNAPSHOT ' num2str(cont) ' value ' num2str([v1(i),v2(j)]) '-------'])
%         data.ip.direction = v2(j)*180/pi;
%         data = run_domain(data,[]);
%         data = run_boundary(data,[]);
%         ufem(:,i,j) = data.solution;
        ufem(:,i,j) = computeFEMfromPGDdata([v1(i),v2(j)],PGDmeshes,parameters);
        cont = cont + 1;
    end
end

bigsave('ufemdataMataro_10x10','ufem')