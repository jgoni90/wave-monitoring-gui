function f = butterfly_A(meshes,mode,matrices)

%Coefs from vectors int{b}, int{c},...
bcoef = mode{2}' * matrices(3).f;
ccoef = mode{3}' * matrices(4).f;
dcoef = mode{4}' * matrices(5).f;
ecoef = mode{5}' * matrices(6).f;

%Coefs from vectors int{Bb}, int{exp(cos(T))}
Bbcoef = meshes(3).X' * matrices(3).M * mode{2};
Ttcoef = mode{1}' * matrices(1).M * exp(cos(meshes(1).X));

%Coefs from non-linear vectors int{cos(CT)tc}, int{sin^D(T/E)tde}
f1 = nonlin1([meshes(1),meshes(4)],mode(3));
nlcoef1 = mode{1}' * f1;
f2 = nonlin2_T([meshes(1),meshes(5),meshes(6)],[mode(4),mode(5)]);
nlcoef2 = mode{1}' * f2;

%Vector f of the current mode A
f = (matrices(2).M * meshes(2).X) * bcoef * Ttcoef * ccoef * dcoef * ecoef -...
    matrices(2).f * nlcoef1 * Bbcoef * dcoef * ecoef +...
    matrices(2).f * nlcoef2 * bcoef * ccoef;