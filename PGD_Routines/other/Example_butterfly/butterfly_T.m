function f = butterfly_T(meshes,mode,matrices)

%Coefs from vectors int{a}, int{b},...
acoef = mode{1}' * matrices(2).f;
bcoef = mode{2}' * matrices(3).f;
ccoef = mode{3}' * matrices(4).f;
dcoef = mode{4}' * matrices(5).f;
ecoef = mode{5}' * matrices(6).f;

%Coefs from vectors int{Aa}, int{Bb}
Aacoef = meshes(2).X' * matrices(2).M * mode{1};
Bbcoef = meshes(3).X' * matrices(3).M * mode{2};

%Coefs from non-linear vectors int{cos(CT)tc}, int{sin^D(T/E)tde}
f1 = nonlin1([meshes(1),meshes(4)],mode(3));
f2 = nonlin2_T([meshes(1),meshes(5),meshes(6)],[mode(4),mode(5)]);

%Vector f of the current mode T
f = (matrices(1).M * exp(cos(meshes(1).X))) * Aacoef * bcoef * ccoef * dcoef * ecoef -...
    f1 * Bbcoef * acoef * dcoef * ecoef +...
    f2 * bcoef * acoef * ccoef;