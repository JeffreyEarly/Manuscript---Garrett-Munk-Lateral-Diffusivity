% plot several examples of signal (u), with vertical axis set to std.
% deviation.

%fraction of HKE from variance to HKE

scaleFactor = 1;
LoadFigureDefaults;

runtype = 'nonlinear';
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
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart_decomp.nc'); 
    matfile = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart_decomp.mat'); 
else
    error('invalid run type.');
end

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

k = ncread(file, 'k');
l = ncread(file, 'l');
j = ncread(file, 'j');

[K,L,J] = ndgrid(k,l,j);
Kh = sqrt(K.*K + L.*L);
RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

nK = length(kAxis)-1;
nModes = length(j);

t = ncread(file, 't');
nT = length(t);

variables = {'Ap_realp', 'Ap_imagp', 'Am_realp', 'Am_imagp', 'B_realp', 'B_imagp'};
Ppm_HKE_factor = wavemodel.Ppm_HKE_factor;
P0_HKE_factor = wavemodel.P0_HKE_factor;
conversion_factor = {Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,Ppm_HKE_factor,P0_HKE_factor,P0_HKE_factor};
Nvars = length(variables);

ncid = netcdf.open(file);
variableIDs = zeros(length(variables),1);
for i=1:length(variables)
    variableIDs(i) = netcdf.inqVarID(ncid,variables{i});
end

% Fully linear wave
iVars = 1:4;
iMode = 8;
iK = 8;
type = 'Wave';
titlestring = sprintf('Linear wave coefficients (k,j) = (%d km,%d)',round(2*pi/kAxis(iK)/1e3),iMode);

% Nonlinear wave
iVars = 1:4;
iMode = 5;
iK = 55;
type = 'Wave';
titlestring = sprintf('Linear wave coefficients (k,j) = (%d km,%d)',round(2*pi/kAxis(iK)/1e3),iMode);

% Vortex in nonlinear wave region
iVars = 5:6;
iMode = 5;
iK = 55;
type = 'Vortex';
titlestring = sprintf('Linear vortex coefficients (k,j) = (%d km,%d)',round(2*pi/kAxis(iK)/1e3),iMode);



ExampleSeries = zeros(length(t),5);

indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(RedundantCoefficients(:,:,1)) );
        
AC = zeros(length(t),1);
DOF = zeros(length(t),1);
nloops = 0;
HKE = 0;
HKEhist = [];
HKEFromVariance = 0;
for iIndex = 1:length(indicesForK)
    [i,j] = ind2sub([size(K,1) size(K,2)],indicesForK(iIndex));
    
    % Need to fix the inertial modes. Somehow I'm not dealing with
    % phases correctly, or something.
    
    for iVar = iVars
        %                 u = double(squeeze(ncread(file, variables{iVar}, [i j iMode 1], [1 1 1 length(t)], [1 1 1 1])));
        u = double(squeeze(netcdf.getVar(ncid, variableIDs(iVar), [i j iMode 1]-1, [1 1 1 length(t)], [1 1 1 1])));
        u = u*sqrt(conversion_factor{iVar}(i,j,iMode));
        [ACu, DOFu] = Autocorrelation(u, length(t)-1);
        
        if any(isnan(ACu))
            continue; % this will occur for the occasional unresolved mode. Seems to only be the Nyquist, which is okay.
        end
        AC = AC + mean(u.*conj(u))*ACu;
        DOF = DOF + DOFu;
        nloops = nloops+1;
        if (nloops<=5)
           ExampleSeries(:,nloops) = u; 
        end
        
        
        
        HKEhist = [HKEhist; mean(u.*conj(u))];
        HKE = HKE + mean(u.*conj(u)); % time mean of the *total* variance of u
        HKEFromVariance = HKEFromVariance + var(u,1); % time variance without time-mean of u
    end
    
end
AC = AC/HKE;
SE =  circshift(sqrt((1 + 2*cumsum(AC.^2,1))./DOF),1);
SE(2) = sqrt(1./DOF(1));
SE(1) = 0;

FigureSize = [50 50 figure_width_2col+8 600*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

subplot(2,1,1)
plot(t/86400,ExampleSeries,'LineWidth', 2), hold on
plot(t/86400,ones(size(t))*2*sqrt(mean(HKEhist)),'Color',0.5*[1 1 1],'LineWidth', 2)
plot(t/86400,-ones(size(t))*2*sqrt(mean(HKEhist)),'Color',0.5*[1 1 1],'LineWidth', 2)
ylim([-1.1 1.1]*max(max(max(ExampleSeries)),2*sqrt(mean(HKEhist))))
ylabel('m^{3/2}/s', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title(titlestring)

subplot(2,1,2)
plot(t/86400,AC,'LineWidth', 2), hold on
plot(t/86400,SE,'k', 'LineWidth', 2)
ylim([0.0 1.1])
ylabel('Correlation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title('Autocorrelation')

% print('-depsc',sprintf('../figures/Autocorrelation%s_k_%d_j_%d.eps',type,iK,iMode))