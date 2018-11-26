% GM
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-26T091749_256x64x65.nc';

% GM with limited K
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-26T111338_256x64x65.nc';

% Smaller GM with limited K, z = -325
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-26T130953_128x32x33.nc';
% Smaller GM with limited K, z = 0
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-26T135610_128x32x33.nc';
% Smaller GM with limited K, z = -350
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-26T142335_128x32x33.nc';

t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nFloatLevels = ncreadatt(file, '/', 'nFloatLevels');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

x = ncread(file, 'x-position-exact')';
y = ncread(file, 'y-position-exact')';
z = ncread(file, 'z-position-exact')';
rho = ncread(file, 'density-exact')';

xLinear = ncread(file, 'x-position-linear')';
yLinear = ncread(file, 'y-position-linear')';
zLinear = ncread(file, 'z-position-linear')';
rhoLinear = ncread(file, 'density-linear')';

xSpline = ncread(file, 'x-position-spline')';
ySpline = ncread(file, 'y-position-spline')';
zSpline = ncread(file, 'z-position-spline')';
rhoSpline = ncread(file, 'density-spline')';

D2_linear = (xLinear-x).^2 + (yLinear-y).^2;
D2_spline = (xSpline-x).^2 + (ySpline-y).^2;

D2z_linear = (zLinear-z).^2;
D2z_spline = (zSpline-z).^2;

kappa_h_linear = D2_linear(end)/(t(end)-t(1))/4;
kappa_h_spline = D2_spline(end)/(t(end)-t(1))/4;

timescale = 1/60;
figure('Name','Horizontal diffusivity')
subplot(2,1,1)
plot(timescale*t(2:end),(1/4)*[D2_linear(2:end)./t(2:end), D2_spline(2:end)./t(2:end)])
legend('linear interp','spline interp')
ylabel('horizontal diffusivity (m^2/s)')
subplot(2,1,2)
plot(timescale*t(2:end), (1/2)*[D2z_linear(2:end)./t(2:end), D2z_spline(2:end)./t(2:end)])
xlabel('minutes')
ylabel('vertical diffusivity (m^2/s)')

return

kappa_z_exact =  ((dz_drho*(rho(2:end)-rho(1))).^2)./t(2:end);
kappa_z_linear =  ((dz_drho*(rhoLinear(2:end)-rhoLinear(1))).^2)./t(2:end);
kappa_z_spline =  ((dz_drho*(rhoSpline(2:end)-rhoSpline(1))).^2)./t(2:end);

figure('Name','Diapycnal diffusivity')
plot(t(2:end),kappa_z_exact), ylog, hold on
plot(t(2:end),kappa_z_linear)
plot(t(2:end),kappa_z_spline)
legend('exact', 'linear', 'spline')
ylabel('diffusivity (m^2/s)')
xlabel('time (s)')
title('Numerical diffusivity in the vertical')