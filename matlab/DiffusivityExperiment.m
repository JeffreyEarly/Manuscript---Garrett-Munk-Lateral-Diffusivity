%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusivityExperiment
%
% This script uses the InternalWaveModel to generate the time-evolution of
% a Garrett-Munk spectrum of internal waves and save the output to a NetCDF
% file. Advects particles as well.
%
% This has been modified to start from Winter's model initial conditions
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% January 21st, 2020   Version 1.0

outputInterval = 15*60;
maxTime = 10.0*86400; %10*outputInterval;
interpolationMethod = 'spline';
shouldOutputEulerianFields = 0;
nLevels = 5;

outputfolder = '/Volumes/Samsung_T5/nsf_iwv';

precision = 'single';

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

wintersModelFile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart';
WM = WintersModel(wintersModelFile);
wavemodel = WM.wavemodel;

[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);

[t_float,x_float,y_float,z_float] = WM.ParticleTrajectories();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Zero out some wavenumbers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_max = 0.0016;
j_max = 47;

wavemodel.Amp_plus(wavemodel.Kh > k_max | wavemodel.J > j_max ) = 0;
wavemodel.Amp_minus(wavemodel.Kh > k_max | wavemodel.J > j_max) = 0;
wavemodel.B(wavemodel.Kh > k_max | wavemodel.J > j_max) = 0;

wavemodel.GenerateWavePhases(wavemodel.Amp_plus,wavemodel.Amp_minus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Grab the floats
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab the intial condition, and convert to column vectors
x_float = reshape(x_float(1,:),[],1);
y_float = reshape(y_float(1,:),[],1);
z_float = reshape(z_float(1,:),[],1)-wavemodel.Lz;

% this the the flux, in dy/dt = f(y,t)
f = @(t,y) wavemodel.VelocityFieldAtTimePosition(t,y,interpolationMethod);

% intial conditions.
xyz0 = cat(2,x_float,y_float,z_float);
nFloats = length(x_float);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the proper time interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

period = 2*pi/wavemodel.N0;
[u,v] = wavemodel.VelocityFieldAtTime(0.0);
U = max(max(max( sqrt(u.*u + v.*v) )));

cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
oscillatoryDT = period/12;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

deltaT = outputInterval/ceil(outputInterval/deltaT);
fprintf('Rounding to match the output interval dt: %.2f\n',deltaT);

t = (0:outputInterval:maxTime)';
if t(end) < period
    t(end+1) = period;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = sprintf('%s/DiffusivityExperiment_%s_%dx%dx%d.nc', outputfolder,datestr(datetime('now'),'yyyy-mm-ddTHHMMSS'),wavemodel.Nx,wavemodel.Ny,wavemodel.Nz);

% Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
totalFields = 4;
totalSize = totalFields*bytePerFloat*length(t)*(wavemodel.Nx)*(wavemodel.Ny)*(wavemodel.Nz)/1e9;
fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',filepath,totalSize);

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

% Define the dynamical variables
if shouldOutputEulerianFields == 1
    uVarID = netcdf.defVar(ncid, 'u', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    vVarID = netcdf.defVar(ncid, 'v', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    wVarID = netcdf.defVar(ncid, 'w', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    zetaVarID = netcdf.defVar(ncid, 'zeta', ncPrecision, [xDimID,yDimID,zDimID,tDimID]);
    netcdf.putAtt(ncid,uVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,vVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,wVarID, 'units', 'm/s');
    netcdf.putAtt(ncid,zetaVarID, 'units', 'm');
end

% Define the *float* dimensions
floatDimID = netcdf.defDim(ncid, 'float_id', nFloats);
xFloatID = netcdf.defVar(ncid, 'x-position', ncPrecision, [floatDimID,tDimID]);
yFloatID = netcdf.defVar(ncid, 'y-position', ncPrecision, [floatDimID,tDimID]);
zFloatID = netcdf.defVar(ncid, 'z-position', ncPrecision, [floatDimID,tDimID]);
densityFloatID = netcdf.defVar(ncid, 'density', ncPrecision, [floatDimID,tDimID]);
netcdf.putAtt(ncid,xFloatID, 'units', 'm');
netcdf.putAtt(ncid,yFloatID, 'units', 'm');
netcdf.putAtt(ncid,zFloatID, 'units', 'm');

% Write some metadata
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', latitude);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', N0);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'GMReferenceLevel', GMReferenceLevel);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from InternalWaveModel.m written by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', wavemodel.version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'max-wavelength-in-spectrum', maxWavelength);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'nFloatLevels', nLevels);
% netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'kappa', kappa);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'interpolation-method', interpolationMethod);

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
integrator = Integrator( f, xyz0, deltaT);
% profile on
for iTime=1:3%length(t)
    if iTime == 2
       startTime = datetime('now'); 
    end
    if iTime == 3 || mod(iTime,10) == 0
        timePerStep = (datetime('now')-startTime)/(iTime-2);
        timeRemaining = (length(t)-iTime+1)*timePerStep;   
        fprintf('\twriting values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iTime, length(t), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end

    if shouldOutputEulerianFields == 1
        [u,v]=wavemodel.VelocityFieldAtTime(t(iTime));
        [w,zeta] = wavemodel.VerticalFieldsAtTime(t(iTime));
    
        netcdf.putVar(ncid, setprecision(uVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], u);
        netcdf.putVar(ncid, setprecision(vVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], v);
        netcdf.putVar(ncid, setprecision(wVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], w);
        netcdf.putVar(ncid, setprecision(zetaVarID), [0 0 0 iTime-1], [wavemodel.Nx wavemodel.Ny wavemodel.Nz 1], zeta);
    end
    netcdf.putVar(ncid, setprecision(tVarID), iTime-1, 1, t(iTime));
    
    p = integrator.StepForwardToTime(t(iTime));
    netcdf.putVar(ncid, setprecision(xFloatID), [0 iTime-1], [nFloats 1], p(:,1));
    netcdf.putVar(ncid, setprecision(yFloatID), [0 iTime-1], [nFloats 1], p(:,2));
    netcdf.putVar(ncid, setprecision(zFloatID), [0 iTime-1], [nFloats 1], p(:,3));
    netcdf.putVar(ncid, setprecision(densityFloatID), [0 iTime-1], [nFloats 1], wavemodel.DensityAtTimePosition(t(iTime),p(:,1),p(:,2),p(:,3),interpolationMethod)-wavemodel.rho0);
end
% profile viewer
fprintf('Ending numerical simulation on %s\n', datestr(datetime('now')));

netcdf.close(ncid);	
