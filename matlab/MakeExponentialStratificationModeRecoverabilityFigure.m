%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% How many modes are resolvable with exponential stratification?
% 
% Notes: for 256 grid points in the vertical you can resolve 67 vertical
% modes. There is no dependence on wavenumber k for the until you get
% around 100 meters.

latitude = 31;

Lx = 800e3;
Lz = 5000;

Nz = 513;
Nx = 128;

dk = 1/Lx;          % fourier frequency
k = 2*pi*([0:ceil(Nx/2)-1 -floor(Nx/2):-1]*dk)';

zIn = [-Lz 0];
z = linspace(-Lz,0,Nz)';

if 1==1
    im = InternalModesExponentialStratification([5.2e-3 1025], zIn, z, latitude,'nModes',floor(Nz/2)); %floor(Nz/3.9)
%     [F,G] = im.ModesAtWavenumber( 0.0565);
else
    N0 = 5.2e-3; % Choose your stratification 7.6001e-04
    
    rho0 = 1025; g = 9.81;
    rho = @(z) -(N0*N0*rho0/g)*z + rho0;
    im = InternalModes(rho,[-Lz 0],z,latitude);
end
% im = InternalModesConstantStratification(5.2e-3, zIn, z, 33);
im.normalization = Normalization.kConstant;
[F,G] = im.ModesAtWavenumber( 0);



nGoodModes_F = InternalModes.NumberOfWellConditionedModes(F);
nGoodModes_G = InternalModes.NumberOfWellConditionedModes(G);

min([nGoodModes_F nGoodModes_G])

kappa_F = InternalModes.ConditionNumberAsFunctionOfModeNumber(F);
kappa_G = InternalModes.ConditionNumberAsFunctionOfModeNumber(G);
figure, plot([kappa_F kappa_G]),ylog, hold on, vlines([nGoodModes_F nGoodModes_G],'g--')