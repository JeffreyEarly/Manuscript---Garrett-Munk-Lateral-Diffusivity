% plot several examples of signal (u), with vertical axis set to std.
% deviation.

%fraction of HKE from variance to HKE

scaleFactor = 1;
LoadFigureDefaults;

NonlinearSpindownFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
NonlinearForcedFromInitialConditionsFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
LinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';
NonlinearSteadyStateFile = '/Volumes/Samsung_T5/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';

dynamicalfile = NonlinearSteadyStateFile;
file = '/Volumes/Samsung_T5/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.nc';
matfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.mat';

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

% Fully linear
iVars = 1:4;
iMode = 2;
iK = 8;
titlestring = sprintf('Linear wave coefficients (k,j) = (%d,%d)',iK,iMode);

% Partiatially nonlinear
iVars = 5:6;
iMode = 27;
iK = 20;
titlestring = sprintf('Linear vortex coefficients (k,j) = (%d,%d)',iK,iMode);

% Fully nonlinear
iVars = 1:4;
iMode = 228;
iK = 12;
titlestring = sprintf('Linear wave coefficients (k,j) = (%d,%d)',iK,iMode);

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
ylim([-1.1 1.1]*max(max(ExampleSeries)))
ylabel('m^3/s^2', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title(titlestring)

subplot(2,1,2)
plot(t/86400,AC,'LineWidth', 2), hold on
plot(t/86400,SE,'k', 'LineWidth', 2)
ylim([0.0 1.1])
ylabel('Correlation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title('Autocorrelation')

print('-depsc',sprintf('Autocorrelation_k_%d_j_%d.eps',iK,iMode))