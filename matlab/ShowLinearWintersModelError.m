runtype = 'linear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end


WM = WintersModel(file);
wavemodel = WM.wavemodel;

nFiles = WM.NumberOf3DOutputFiles;
[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);

error2 = @(u,u_unit) abs((u-u_unit))/max(abs(u_unit(:)));

fileIncrements = 1:50:nFiles;
t_error=zeros(length(fileIncrements),1);
u_error=zeros(length(fileIncrements),1);
v_error=zeros(length(fileIncrements),1);
w_error=zeros(length(fileIncrements),1);
rho_error=zeros(length(fileIncrements),1);

for iTime = 1:length(fileIncrements)
    [t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(fileIncrements(iTime),'t','u','v','w','rho_prime');
    [u_back,v_back,w_back,rho_prime_back] = wavemodel.VariableFieldsAtTime(t,'u','v','w','rho_prime');
    
    t_error(iTime) = t;
    
    error = error2(u,u_back);
    u_error(iTime) = max(error(:));
    
    error = error2(v,v_back);
    v_error(iTime) = max(error(:));
    
    error = error2(w,w_back);
    w_error(iTime) = max(error(:));
    
    error = error2(rho_prime,rho_prime_back);
    rho_error(iTime) = max(error(:));
end

figure
plot(t_error,[u_error, v_error, w_error, rho_error])
legend('u','v','w','rho')

w = wavemodel;
wavemodel2 = InternalWaveModelConstantStratification([w.Lx, w.Ly, w.Lz], [w.Nx, w.Ny, w.Nz], w.latitude, w.N0, w.rho0);
wavemodel2.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);