%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LagrangianErrorExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. It advects particles using the exact (spectral) interpolation,
% linear interpolation and cubic-spline interpolation. The goal is to
% assess the errors of these techniques.
%
% November 28th, 2018 tests.
% Linear interpolation does very well if you zero-pad with half of the
% wavenumbers, but cubic spline does even better that it's probably worth
% the 33% penalty in time, just to have the extra accuracy. Basically I'm
% finding a numerical diffusivity of O(1e-4) using cubic splines and
% O(1e-2) using linear interpolation.
%
% The other big question is what time step can you get away with. The CFL
% condition requires either the advective speed (~20 cm/s) or the wave speed
% (2 m/s). I find that in *low* resolution situations, you should set the
% time step to the cfl=1 for the wave speed. But, as I get to higher
% resolution, 0.5 cfl for the advective time step seems fine. I don't know
% why this is! Probably because particle motion at the smallest scales
% changes what is controlling it?
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% May 2nd, 2017         Version 1.0
% November 28th, 2018   Version 2.0

N = 16;
aspectRatio = 4;
L = 25e3;

Lx = aspectRatio*L;
Ly = L;
Lz = 1300;

Nx = aspectRatio*N;
Ny = N;
Nz = N+1; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification
GMReferenceLevel = 1.0;

outputInterval = 1440;
maxTime = 86400;

outputfolder = '/Volumes/Samsung_T5/linear-model-diffusivity';

precision = 'double';

if strcmp(precision,'single')
    ncPrecision = 'NC_FLOAT';
    setprecision = @(x) single(x);
    bytePerFloat = 4;
else
    ncPrecision = 'NC_DOUBLE';
    setprecision = @(x) double(x);
    bytePerFloat = 8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shouldUseGMSpectrum = 0;

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
maxK = max(0.25*abs(wavemodel.k));
maxMode = max(round(0.25*wavemodel.j));
wavemodel.InitializeWithGMSpectrum(GMReferenceLevel,'maxK',maxK,'maxMode',maxMode);
        
wavemodel.ShowDiagnostics();
period = 2*pi/wavemodel.N0;
[u,v] = wavemodel.VelocityFieldAtTime(0.0);
U = max(max(max( sqrt(u.*u + v.*v) )));
c = max(max(max(wavemodel.C)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create floats/drifters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = wavemodel.x(2)-wavemodel.x(1);
dy = wavemodel.y(2)-wavemodel.y(1);
N = 1;
nLevels = 1;
x_float = (1:N)*dx;
y_float = (1:N)*dy;
z_float = -650;

[x_float,y_float,z_float] = ndgrid(x_float,y_float,z_float);
x_float = reshape(x_float,[],1);
y_float = reshape(y_float,[],1);
z_float = reshape(z_float,[],1);
nFloats = numel(x_float);

% Now let's place the floats along an isopycnal.
isopycnalDeviation = wavemodel.IsopycnalDisplacementAtTimePosition(0,x_float,y_float,z_float, 'exact');
z_isopycnal = z_float + isopycnalDeviation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfl = 0.25;
advectiveDT = (wavemodel.x(2)-wavemodel.x(1))/U;
waveDT = (wavemodel.x(2)-wavemodel.x(1))/c;
if advectiveDT < waveDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the wave dt: %.2f\n',waveDT);
    deltaT = waveDT;
end

deltaT = outputInterval/ceil(outputInterval/deltaT);
fprintf('Rounding to match the output interval dt: %.2f\n',deltaT);
fprintf('advective cfl: %.2f, wave cfl: %.2f\n',deltaT/advectiveDT, deltaT/waveDT);

t_in = (0:outputInterval:maxTime)';
p0 = [x_float, y_float, z_float];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/LagrangianErrorExperiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),Nx,Ny,Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 3;
totalSize = totalFields*bytePerFloat*length(t_in)*3*nFloats/1e6;
fprintf('Writing output file to %s\nExpected file size is %.2f MB.\n',filepath,totalSize);

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
ncid = netcdf.create(filepath, cmode);

% Define the dimensions
xDimID = netcdf.defDim(ncid, 'x', wavemodel.Nx);
yDimID = netcdf.defDim(ncid, 'y', wavemodel.Ny);
zDimID = netcdf.defDim(ncid, 'z', wavemodel.Nz);
tDimID = netcdf.defDim(ncid, 't', netcdf.getConstant('NC_UNLIMITED'));

% Define the coordinate variables
xVarID = netcdf.defVar(ncid, 'x', ncPrecision, xDimID);
yVarID = netcdf.defVar(ncid, 'y', ncPrecision, yDimID);
zVarID = netcdf.defVar(ncid, 'z', ncPrecision, zDimID);
tVarID = netcdf.defVar(ncid, 't', ncPrecision, tDimID);
netcdf.putAtt(ncid,xVarID, 'units', 'm');
netcdf.putAtt(ncid,yVarID, 'units', 'm');
netcdf.putAtt(ncid,zVarID, 'units', 'm');
netcdf.putAtt(ncid,tVarID, 'units', 's');

% Define the *float* dimensions --- spectral interpolation, adaptive RK45
floatDimID = netcdf.defDim(ncid, 'float_id', nFloats);
xFloatID = netcdf.defVar(ncid, 'x-position-exact', ncPrecision, [floatDimID,tDimID]);
yFloatID = netcdf.defVar(ncid, 'y-position-exact', ncPrecision, [floatDimID,tDimID]);
zFloatID = netcdf.defVar(ncid, 'z-position-exact', ncPrecision, [floatDimID,tDimID]);
densityFloatID = netcdf.defVar(ncid, 'density-exact', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFloatID, 'units', 'm');

% Define the *float* dimensions --- spectral interpolation, fixed dt RK4
xFixedDtFloatID = netcdf.defVar(ncid, 'x-position-exact-fixed-dt', ncPrecision, [floatDimID,tDimID]);
yFixedDtFloatID = netcdf.defVar(ncid, 'y-position-exact-fixed-dt', ncPrecision, [floatDimID,tDimID]);
zFixedDtFloatID = netcdf.defVar(ncid, 'z-position-exact-fixed-dt', ncPrecision, [floatDimID,tDimID]);
densityFixedDtFloatID = netcdf.defVar(ncid, 'density-exact-fixed-dt', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFixedDtFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFixedDtFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFixedDtFloatID, 'units', 'm');

% Define the *float* dimensions --- linear interpolation, adaptive RK45
xLinearFloatID = netcdf.defVar(ncid, 'x-position-linear', ncPrecision, [floatDimID,tDimID]);
yLinearFloatID = netcdf.defVar(ncid, 'y-position-linear', ncPrecision, [floatDimID,tDimID]);
zLinearFloatID = netcdf.defVar(ncid, 'z-position-linear', ncPrecision, [floatDimID,tDimID]);
densityLinearFloatID = netcdf.defVar(ncid, 'density-linear', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,yLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,zLinearFloatID, 'units', 'm');

% Define the *float* dimensions --- linear interpolation, fixed dt RK4
xFixedDtLinearFloatID = netcdf.defVar(ncid, 'x-position-linear-fixed-dt', ncPrecision, [floatDimID,tDimID]);
yFixedDtLinearFloatID = netcdf.defVar(ncid, 'y-position-linear-fixed-dt', ncPrecision, [floatDimID,tDimID]);
zFixedDtLinearFloatID = netcdf.defVar(ncid, 'z-position-linear-fixed-dt', ncPrecision, [floatDimID,tDimID]);
densityFixedDtLinearFloatID = netcdf.defVar(ncid, 'density-linear-fixed-dt', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFixedDtLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFixedDtLinearFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFixedDtLinearFloatID, 'units', 'm');

% Define the *float* dimensions --- spline interpolation, adaptive RK45
xSplineFloatID = netcdf.defVar(ncid, 'x-position-spline', ncPrecision, [floatDimID,tDimID]);
ySplineFloatID = netcdf.defVar(ncid, 'y-position-spline', ncPrecision, [floatDimID,tDimID]);
zSplineFloatID = netcdf.defVar(ncid, 'z-position-spline', ncPrecision, [floatDimID,tDimID]);
densitySplineFloatID = netcdf.defVar(ncid, 'density-spline', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xSplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,ySplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,zSplineFloatID, 'units', 'm');

% Define the *float* dimensions --- spline interpolation, fixed dt RK4
xFixedDtSplineFloatID = netcdf.defVar(ncid, 'x-position-spline-fixed-dt', ncPrecision, [floatDimID,tDimID]);
yFixedDtSplineFloatID = netcdf.defVar(ncid, 'y-position-spline-fixed-dt', ncPrecision, [floatDimID,tDimID]);
zFixedDtSplineFloatID = netcdf.defVar(ncid, 'z-position-spline-fixed-dt', ncPrecision, [floatDimID,tDimID]);
densityFixedDtSplineFloatID = netcdf.defVar(ncid, 'density-spline-fixed-dt', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFixedDtSplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFixedDtSplineFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFixedDtSplineFloatID, 'units', 'm');

% Write some metadata
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', N0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'nFloatLevels', nLevels);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'maxK', maxK);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'maxMode', maxMode);

% End definition mode
netcdf.endDef(ncid);

% Add the data for the coordinate variables
netcdf.putVar(ncid, setprecision(xVarID), wavemodel.x);
netcdf.putVar(ncid, setprecision(yVarID), wavemodel.y);
netcdf.putVar(ncid, setprecision(zVarID), wavemodel.z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the model, and write the output to NetCDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startTime = datetime('now');
fprintf('Starting numerical simulation on %s\n', datestr(startTime));

netcdf.putVar(ncid, tVarID, 0, length(t_in), t_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Spectral interpolation with adaptive time step...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'exact');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xFloatID, x);
netcdf.putVar(ncid, yFloatID, y);
netcdf.putVar(ncid, zFloatID, z);
netcdf.putVar(ncid, densityFloatID, setprecision(rho));
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-exact-adaptive', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Spectral interpolation with fixed time step starting on %s...',datestr(startTime));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integrator = Integrator( f, p0, deltaT);
t = t_in;
for iTime=1:length(t)
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, xFixedDtFloatID, [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, yFixedDtFloatID, [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, zFixedDtFloatID, [0 iTime-1], [nFloats 1], p(:,3));
    netcdf.putVar(ncid, densityFixedDtFloatID, [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,1),p(:,2),p(:,3), 'exact')-wavemodel.rho0);
end
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-exact-fixed', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Linear interpolation with adaptive time step starting on %s...', startTime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'linear');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xLinearFloatID, setprecision(x));
netcdf.putVar(ncid, yLinearFloatID, setprecision(y));
netcdf.putVar(ncid, zLinearFloatID, setprecision(z));
netcdf.putVar(ncid, densityLinearFloatID, setprecision(rho));
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-linear-adaptive', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Linear interpolation with fixed time step starting on %s...',datestr(startTime));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integrator = Integrator( f, p0, deltaT);
t = t_in;
for iTime=1:length(t)
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, xFixedDtLinearFloatID, [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, yFixedDtLinearFloatID, [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, zFixedDtLinearFloatID, [0 iTime-1], [nFloats 1], p(:,3));
    netcdf.putVar(ncid, densityFixedDtLinearFloatID, [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,1),p(:,2),p(:,3), 'exact')-wavemodel.rho0);
end
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-linear-fixed', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Spline interpolation with adaptive time step starting on %s...', datestr(startTime));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(t,y) wavemodel.VelocityAtTimePositionVector(t,y,'spline');
[t,p] = ode45(f,t_in, p0,odeset('RelTol',1e-11,'AbsTol',1e-8));
x = p(:,1);
y = p(:,2);
z = p(:,3);
rho = zeros(size(x));
for iTime=1:length(t)
    rho(iTime) = wavemodel.DensityAtTimePosition(t(iTime),x(iTime),y(iTime),z(iTime), 'exact') - wavemodel.rho0;
end
netcdf.putVar(ncid, xSplineFloatID, setprecision(x));
netcdf.putVar(ncid, ySplineFloatID, setprecision(y));
netcdf.putVar(ncid, zSplineFloatID, setprecision(z));
netcdf.putVar(ncid, densitySplineFloatID, setprecision(rho));
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-spline-adaptive', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Spline interpolation with fixed time step starting on %s...',datestr(startTime));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integrator = Integrator( f, p0, deltaT);
t = t_in;
for iTime=1:length(t)
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, xFixedDtSplineFloatID, [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, yFixedDtSplineFloatID, [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, zFixedDtSplineFloatID, [0 iTime-1], [nFloats 1], p(:,3));
    netcdf.putVar(ncid, densityFixedDtSplineFloatID, [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,1),p(:,2),p(:,3), 'exact')-wavemodel.rho0);
end
integrationTime = round(seconds(datetime('now')-startTime));
netcdf.reDef(ncid)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'integration-time-spline-fixed', integrationTime);
netcdf.endDef(ncid);
fprintf('total integration time %d seconds.\n',integrationTime);

netcdf.close(ncid);	
