% GM
% high resolution, nearly full spectrum model run that includes fixed time
% step advected particles
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T063305_256x64x65.nc';

% trying to take too big of a time step to watch it fall apart
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T094611_64x16x17.nc';

% trying to take too big of a time step to watch it fall apart
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T094928_64x16x17.nc';

% trying to take too big of a time step to watch it fall apart
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T095147_64x16x17.nc';

% trying to take too big of a time step to watch it fall apart -- advective
% cfl = 1/2 and finally things fell apart
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T095856_64x16x17.nc';

% trying to take too big of a time step to watch it fall apart -- advective
% cfl back down to 1/4 still badd
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T100510_64x16x17.nc';

% trying to take too big of a time step to watch it fall apart -- advective
% cfl back down to 1/8 -- but wave cfl is 0.8
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T100950_64x16x17.nc';

% advective cfl: 0.09, wave cfl: 0.99 --- all works well. Let's try to
% break this again.
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T101349_128x32x33.nc';

% advective cfl: advective cfl: 0.18, wave cfl: 1.98 --- starting to show
% signs of weakness, but still did okay.
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T103141_128x32x33.nc';

% advective cfl: advective cfl: 0.44, wave cfl: 3.97 --- yes, fell apart,
% but just barely.
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T104838_128x32x33.nc';

% advective cfl: advective cfl: 0.45, wave cfl: 3.97 --- bumping
% resolution. Can we pull back together? Yes. Things are stable again.
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T113243_256x64x65.nc';

% now using 1/4 of the wavenumbers
file = '/Volumes/Samsung_T5/linear-model-diffusivity/LagrangianErrorExperiment_2018-11-28T150337_64x16x17.nc';


t = ncread(file, 't');

Nx = length(ncread(file, 'x'));
Ny = length(ncread(file, 'y'));
Nz = length(ncread(file, 'z'));

nFloatLevels = ncreadatt(file, '/', 'nFloatLevels');
N0 = ncreadatt(file, '/', 'N0');
rho0 = 1025;
dz_drho = 9.81/(N0*N0*rho0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

xFixedDt = ncread(file, 'x-position-exact-fixed-dt')';
yFixedDt = ncread(file, 'y-position-exact-fixed-dt')';
zFixedDt = ncread(file, 'z-position-exact-fixed-dt')';
rhoFixedDt = ncread(file, 'density-exact-fixed-dt')';

xFixedDtLinear = ncread(file, 'x-position-linear-fixed-dt')';
yFixedDtLinear = ncread(file, 'y-position-linear-fixed-dt')';
zFixedDtLinear = ncread(file, 'z-position-linear-fixed-dt')';
rhoFixedDtLinear = ncread(file, 'density-linear-fixed-dt')';

xFixedDtSpline = ncread(file, 'x-position-spline-fixed-dt')';
yFixedDtSpline = ncread(file, 'y-position-spline-fixed-dt')';
zFixedDtSpline = ncread(file, 'z-position-spline-fixed-dt')';
rhoFixedDtSpline = ncread(file, 'density-spline-fixed-dt')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D2_linear = (xLinear-x).^2 + (yLinear-y).^2;
D2_spline = (xSpline-x).^2 + (ySpline-y).^2;
D2_exact_fixed_dt = (xFixedDt-x).^2 + (yFixedDt-y).^2;
D2_linear_fixed_dt = (xFixedDtLinear-x).^2 + (yFixedDtLinear-y).^2;
D2_spline_fixed_dt = (xFixedDtSpline-x).^2 + (yFixedDtSpline-y).^2;

D2z_linear = (zLinear-z).^2;
D2z_spline = (zSpline-z).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa_h_linear = D2_linear(2:end)./t(2:end)/4;
kappa_h_spline = D2_spline(2:end)./t(2:end)/4;
kappa_h_exact_fixed_dt = D2_exact_fixed_dt(2:end)./t(2:end)/4;
kappa_h_linear_fixed_dt = D2_linear_fixed_dt(2:end)./t(2:end)/4;
kappa_h_spline_fixed_dt = D2_spline_fixed_dt(2:end)./t(2:end)/4;

timescale = 1/60;
figure('Name','Horizontal diffusivity')
subplot(2,1,1)
plot(timescale*t(2:end),[kappa_h_linear,kappa_h_spline,kappa_h_exact_fixed_dt,kappa_h_linear_fixed_dt,kappa_h_spline_fixed_dt])
legend('linear interp','spline interp','spectral-fixed-dt','linear-fixed-dt','spline-fixed dt')
ylabel('horizontal diffusivity (m^2/s)'), ylog
subplot(2,1,2)
plot(timescale*t(2:end), (1/2)*[D2z_linear(2:end)./t(2:end), D2z_spline(2:end)./t(2:end)])
ylog
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