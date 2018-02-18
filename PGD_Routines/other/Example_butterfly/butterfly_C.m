function f = butterfly_C(meshes,mode,matrices)

%Coefs from vectors int{a}, int{b},...
acoef = mode{2}' * matrices(2).f;
bcoef = mode{3}' * matrices(3).f;
dcoef = mode{4}' * matrices(5).f;
ecoef = mode{5}' * matrices(6).f;

%Coefs from vectors int{Aa}, int{Bb}, int{exp(cos(T))}
Aacoef = meshes(2).X' * matrices(2).M * mode{2};
Bbcoef = meshes(3).X' * matrices(3).M * mode{3};
Ttcoef = mode{1}' * matrices(1).M * exp(cos(meshes(1).X));

%Coefs from non-linear vectors int{cos(CT)tc}, int{sin^D(T/E)tde}
f1 = nonlin1([meshes(4),meshes(1)],mode(1));
f2 = nonlin2_T([meshes(1),meshes(5),meshes(6)],[mode(4),mode(5)]);
nlcoef2 = mode{1}' * f2;

%Vector f of the current mode C
f = matrices(4).f * Aacoef * bcoef * Ttcoef * dcoef * ecoef -...
    f1 * Bbcoef * acoef * dcoef * ecoef +...
    matrices(4).f * nlcoef2 * bcoef * acoef;