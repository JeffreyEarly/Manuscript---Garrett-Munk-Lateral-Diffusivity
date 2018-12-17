runtype = 'linear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
end

if strcmp(runtype,'linear')
    dynamicalfile = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    dynamicalfile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end
output_directory = baseURL;

[filepath,name,ext] = fileparts(dynamicalfile);
file = fullfile(output_directory,strcat(name,'_decomp.nc'));

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

k = ncread(file, 'k');
l = ncread(file, 'l');
j = ncread(file, 'j');

t = ncread(file, 't');
nT = length(t);

variables = {'Ap_realp', 'Ap_imagp', 'Am_realp', 'Am_imagp', 'B_realp', 'B_imagp'};

ncid = netcdf.open(file);
variableIDs = zeros(length(variables),1);
for iK=1:length(variables)
    variableIDs(iK) = netcdf.inqVarID(ncid,variables{iK});
end

iVar = 1; iK = 0; iL=0; iMode = 0;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
Ap_inertial = u+sqrt(-1)*v;

iVar = 3;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
Am_inertial = u+sqrt(-1)*v;

iVar = 5;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
B_inertial = u+sqrt(-1)*v;

A_inertial = Ap_inertial + Am_inertial;

figure('Name',sprintf('Amp and phase of inertial mode'))
subplot(2,1,1)
plot(t,abs(Ap_inertial))
title(sprintf('Amp and phase of inertial mode'))
subplot(2,1,2)
plot(t,angle(Ap_inertial))


iVar = 1; iK = 5; iL=0; iMode = 0;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
Ap_tidal= u+sqrt(-1)*v;

iVar = 3;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
Am_tidal = u+sqrt(-1)*v;

iVar = 5;
u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [iK iL iMode 0], [1 1 1 length(t)])));
v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [iK iL iMode 0], [1 1 1 length(t)])));
B_tidal = u+sqrt(-1)*v;

figure('Name',sprintf('Amp and phase of mode (k,l,j)=(%d,%d,%d)',iK,iL,iMode+1))
subplot(2,1,1)
plot(t,abs(Am_tidal)), hold on, plot(t,abs(Ap_tidal))
title(sprintf('Amp and phase of mode (k,l,j)=(%d,%d,%d)',iK,iL,iMode+1))
subplot(2,1,2)
plot(t,angle(Am_tidal)), hold on, plot(t,angle(Ap_tidal))

return




for iVar = 1:2:length(variableIDs)
    
end