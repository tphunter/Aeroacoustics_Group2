%% Tonal Noise Functions

sres = param.sres;
phi = param.phi
for m = 1:1:20
    for s = -sres:1:sres
        Pmb = i*m*B^2*Omega/(4*pi*c*R0)*