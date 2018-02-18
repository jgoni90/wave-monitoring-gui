function f = butterfly_E(meshes,mode,matrices)

%Coefs from vectors int{a}, int{b},...
acoef = mode{2}' * matrices(2).f;
bcoef = mode{3}' * matrices(3).f;
ccoef = mode{4}' * matrices(4).f;
dcoef = mode{5}' * matrices(5).f;

%Coefs from vectors int{Aa}, int{Bb}, int{exp(cos(T))}
Aacoef = meshes(2).X' * matrices(2).M * mode{2};
Bbcoef = meshes(3).X' * matrices(3).M * mode{3};
Ttcoef = mode{1}' * matrices(1).M * exp(cos(meshes(1).X));

%Coefs from non-linear vectors int{cos(CT)tc}, int{sin^D(T/E)tde}
f1 = nonlin1([meshes(1),meshes(4)],mode(4));
nlcoef1 = mode{1}' * f1;
f2 = nonlin2_E([meshes(6),meshes(1),meshes(5)],[mode(1),mode(5)]);

%Vector f of the current mode E
f = matrices(6).f * Aacoef * Ttcoef * bcoef * ccoef * dcoef -...
    matrices(6).f * Bbcoef * acoef * nlcoef1 * dcoef +...
    f2 * bcoef * acoef * ccoef;