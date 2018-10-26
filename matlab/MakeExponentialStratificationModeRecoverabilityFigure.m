%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% How many modes are resolvable with exponential stratification?
% 
% Notes: for 256 grid points in the vertical you can resolve 67 vertical
% modes. There is no dependence on wavenumber k for the until you get
% around 100 meters.

Lx = 800e3;
Lz = 4000;

Nz = 256;
Nx = 128;

dk = 1/Lx;          % fourier frequency
k = 2*pi*([0:ceil(Nx/2)-1 -floor(Nx/2):-1]*dk)';

zIn = [-Lz 0];
z = linspace(-Lz,0,Nz)';

im = InternalModesExponentialStratification([5.2e-3 1025], zIn, z, 33,'nModes',max(Nz/2,128));
im.normalization = Normalization.kConstant;
[F,G] = im.ModesAtWavenumber(1e-2);

kappa_F = InternalModes.ConditionNumberAsFunctionOfModeNumber(F);
kappa_G = InternalModes.ConditionNumberAsFunctionOfModeNumber(G);

nGoodModes_F = InternalModes.NumberOfWellConditionedModes(F);
nGoodModes_G = InternalModes.NumberOfWellConditionedModes(G);

min([nGoodModes_F nGoodModes_G])

figure, plot([kappa_F kappa_G]),ylog, hold on, vlines([nGoodModes_F nGoodModes_G],'g--')