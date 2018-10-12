NonlinearSpindownFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
NonlinearForcedFromInitialConditionsFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
LinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';
NonlinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';

dynamicalfile = NonlinearSteadyStateFile;
WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

[u,v] = WM.VariableFieldsFrom3DOutputFileAtIndex(10,'u','v');

Nx = wavemodel.Nx;
Ny = wavemodel.Ny;
Nz = wavemodel.Nz;

% Compute the power spectra
u_bar = fft(fft(u,Nx,1),Ny,2)/Nx/Ny;
S_u = u_bar .* conj(u_bar) * wavemodel.Lx * wavemodel.Ly; % m^2 * m^2 /s^2
v_bar = fft(fft(v,Nx,1),Ny,2)/Nx/Ny;
S_v = v_bar .* conj(v_bar) * wavemodel.Lx * wavemodel.Ly; % m^2 * m^2 /s^2


deltaK = ones(size(u_bar)) * (wavemodel.k(2)-wavemodel.k(1));
S_u = S_u.*deltaK;
S_v = S_v.*deltaK;

Kh = wavemodel.Kh;
dK = wavemodel.k(2)-wavemodel.k(1);
kMag = 0:dK:max(max(max(Kh)));
indices = cell( length(kMag), 1);
for i = 1:length(kMag)
    indices{i} = find( squeeze(Kh(:,:,1)) >= kMag(i)-dK/2 & squeeze(Kh(:,:,1)) < kMag(i)+dK/2 );
end

S_u_isotropic = zeros(length(kMag),Nz);
S_v_isotropic = zeros(length(kMag),Nz);
for iK = 1:length(kMag)
    for iJ = 1:Nz
        ulvl = squeeze(S_u(:,:,iJ));
        S_u_isotropic(iK,iJ) = sum(ulvl(indices{iK}));
        
        vlvl = squeeze(S_v(:,:,iJ));
        S_v_isotropic(iK,iJ) = sum(vlvl(indices{iK}));
    end
end

depth = 1;

figure
plot(kAxis,S_u_isotropic(:,depth)), xlog, ylog
hold on, plot(kAxis,S_u_isotropic(:,128))

iStart = floor(length(kAxis)/2);
iIntersect = floor(length(kAxis)*3/4);
x = kAxis(iIntersect+10);
y = S_u_isotropic(iIntersect,depth);
% y = A x^{-m}, log(y) = log(A) - m*log(x)
A = y*x^12;
damping = A*kAxis.^(-12);
index = find(damping<max(S_u_isotropic(:,depth))/10,1,'first');

hold on, plot(kAxis(index:end),damping(index:end),'k','LineWidth',2)
vlines(max(wavemodel.k),'k--')

% Now let's compute how much amplitude is lost with the given damping after
% the entire run, as a function of wavenumber.

T = 36000;
p = 3;
dx = wavemodel.x(2)-wavemodel.x(1);
nu = power(dx/pi,2*p) / T;
T_total = 10*86400;
dissipated = exp( -2*nu*power(kAxis,2*p)*T_total);
figure, plot(kAxis,1-dissipated), xlog, ylog

