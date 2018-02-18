function f = butterfly_B(meshes,mode,matrices)

%Coefs from vectors int{a}, int{c},...
acoef = mode{2}' * matrices(2).f;
ccoef = mode{3}' * matrices(4).f;
dcoef = mode{4}' * matrices(5).f;
ecoef = mode{5}' * matrices(6).f;

%Coefs from vectors int{Aa}, int{exp(cos(T))}
Aacoef = meshes(2).X' * matrices(2).M * mode{2};
Ttcoef = mode{1}' * matrices(1).M * exp(cos(meshes(1).X));

%Coefs from non-linear vectors int{cos(CT)tc}, int{sin^D(T/E)tde}
f1 = nonlin1([meshes(1),meshes(4)],mode(3));
nlcoef1 = mode{1}' * f1;
f2 = nonlin2_T([meshes(1),meshes(5),meshes(6)],[mode(4),mode(5)]);
nlcoef2 = mode{1}' * f2;

%Vector f of the current mode B
f = matrices(3).f * Aacoef * Ttcoef * ccoef * dcoef * ecoef -...
    (matrices(3).M * meshes(3).X) * nlcoef1 * acoef * dcoef * ecoef +...
    matrices(3).f * nlcoef2 * acoef * ccoef;