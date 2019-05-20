runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
    baseURL = '/Users/jearly/Documents/nsf_iw/';
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2019_05/';
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

outputfile = strcat(baseURL,'WaveVortexEnergyVsTimeSpinup.mat');

if ~exist(outputfile,'file')
    k = ncread(file, 'k');
    l = ncread(file, 'l');
    j = ncread(file, 'j');
    
    t = ncread(file, 't');
    nT = length(t);
    
    variables = {'Ap_realp', 'Ap_imagp', 'Am_realp', 'Am_imagp', 'B_realp', 'B_imagp', 'B0_realp', 'B0_imagp'};
    
    ncid = netcdf.open(file);
    variableIDs = zeros(length(variables),1);
    for iK=1:length(variables)
        variableIDs(iK) = netcdf.inqVarID(ncid,variables{iK});
    end
    
    VortexEnergyConversionFactor = wavemodel.B_HKE_factor + wavemodel.B_PE_factor;
    Vortex0EnergyConversionFactor = wavemodel.B0_HKE_factor;
    
    WaveEnergy_plus = zeros(nT,1);
    WaveEnergy_minus = zeros(nT,1);
    VortexEnergy = zeros(nT,1);
    Vortex0Energy = zeros(nT,1);
    
    for iTime=1:nT
        fprintf('%d..',iTime)
        iVar=1;
        u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        A2 = u.*u + v.*v;
        WaveEnergy_plus(iTime) = sum(sum(sum(A2(:))));
        
        iVar=3;
        u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        A2 = u.*u + v.*v;
        WaveEnergy_minus(iTime) = sum(sum(sum(A2(:))));
        
        iVar=5;
        u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [0 0 0 iTime-1], [length(k) length(l) length(j) 1])));
        B2 = VortexEnergyConversionFactor.*(u.*u + v.*v);
        VortexEnergy(iTime) = sum(sum(sum(B2(:))));
        
        iVar=7;
        u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [0 0 iTime-1], [length(k) length(l) 1])));
        v = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar+1), [0 0 iTime-1], [length(k) length(l) 1])));
        B2 = Vortex0EnergyConversionFactor.*(u.*u + v.*v);
        Vortex0Energy(iTime) = sum(sum(sum(B2(:))));
    end
    
%     save(outputfile,'t','WaveEnergy_plus','WaveEnergy_minus','Vortex0Energy','VortexEnergy')
end

figure
subplot(2,1,1)
plot(t/86400,WaveEnergy_plus), hold on
plot(t/86400,WaveEnergy_minus)
legend('\omega_+','\omega_-')
ylabel('energy (m^3/s^2)')
subplot(2,1,2)
plot(t/86400,VortexEnergy), hold on
plot(t/86400,Vortex0Energy)
legend('vortex', 'mean')
xlabel('days')
ylabel('energy (m^3/s^2)')
packfig(2,1)

print('-depsc','../data/2019_05/figures_EarlyV2_GM_NL_forced_damped_restart/WaveVortexEnergyVsTime.eps')
