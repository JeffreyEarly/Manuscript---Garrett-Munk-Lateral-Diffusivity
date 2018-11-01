runtype = 'linear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_11/';
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

s2 = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'s2');

dx = WM.wavemodel.x(2)-WM.wavemodel.x(1);
dy = WM.wavemodel.y(2)-WM.wavemodel.y(1);
dz = WM.wavemodel.z(2)-WM.wavemodel.z(1);
X = WM.wavemodel.X;
Y = WM.wavemodel.Y;
Z = WM.wavemodel.Z;

dV = dx*dy*dz;



nFiles = WM.NumberOf3DOutputFiles;
fileIncrements = 1:1:nFiles;

nT = length(fileIncrements);
t = zeros(nT,1);
m = zeros(nT,1);
m_x = zeros(nT,1);
m_y = zeros(nT,1);
m_z = zeros(nT,1);
m_xx = zeros(nT,1);
m_yy = zeros(nT,1);
m_zz = zeros(nT,1);
startTime = datetime('now');
for iFile = 1:length(fileIncrements)
    if mod(iFile,10) == 2
        timePerStep = (datetime('now')-startTime)/(iFile-1);
        timeRemaining = (nT-iFile+1)*timePerStep;
        fprintf('reading values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iFile, nT, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end
    
    file = WM.PathOf3DOutputFileAtIndex(fileIncrements(iFile));
    s2 = WM.VariableFieldsFrom3DOutputFileAtIndex(iFile,'s2');
    
    t(iFile) = double(ncread(file,'time'));
    
    m(iFile) = sum(sum(sum(s2)))*dV;
    m_x(iFile) = sum(sum(sum(X.*s2)))*dV/m(iFile);
    m_y(iFile) = sum(sum(sum(Y.*s2)))*dV/m(iFile);
    m_z(iFile) = sum(sum(sum(Z.*s2)))*dV/m(iFile);
    
    Xc = X-m_x(iFile);
    Yc = Y-m_y(iFile);
    Zc = Z-m_z(iFile);
    m_xx(iFile) = sum(sum(sum(Xc.*Xc.*s2)))*dV/m(iFile);
    m_yy(iFile) = sum(sum(sum(Yc.*Yc.*s2)))*dV/m(iFile);
    m_zz(iFile) = sum(sum(sum(Zc.*Zc.*s2)))*dV/m(iFile);
end

D2 = (m_xx+m_yy)/2;
[p,~,mu]=polyfit(t,D2,1);
kappa_xy = (p(1)/mu(2))/2;

[p,~,mu]=polyfit(t,m_zz,1);
kappa_z = (p(1)/mu(2))/2;
            
figure('Name','LateralDiffusivityOfTracer')
subplot(3,1,1)
plot(t/86400,m/m(1))
ylabel('total mass')
subplot(3,1,2)
plot(t/86400,[m_x, m_y]/1e3)
legend('x','y')
ylabel('km')
subplot(3,1,3)
plot(t/86400,[m_xx, m_yy,D2]/1e6)
legend('xx','yy','(xx+yy)/2')
title(sprintf('isotropic diffusivity %.2f m^2/s at scale %.0f km',kappa_xy, sqrt(D2(1))/1e3 ));
ylabel('km^2')
print('-depsc',sprintf('TracerLateral-%s.eps',runtype))

figure('Name','VerticalDiffusivityOfTracer')
subplot(2,1,1)
plot(t/86400,m_z)
ylabel('m')
subplot(2,1,2)
plot(t/86400,m_zz)
title(sprintf('isotropic diffusivity %.2g m^2/s at scale %.2f m',kappa_z, sqrt(m_zz(1)) ));
ylabel('m^2')
print('-depsc',sprintf('TracerVertical-%s.eps',runtype))