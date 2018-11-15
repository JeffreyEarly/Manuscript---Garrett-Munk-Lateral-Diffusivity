scaleFactor = 1;
LoadFigureDefaults;

mooringsfile = '../data/2018_11/EarlyV2_GM_NL_forced_damped_restart_moorings.mat';
% mooringsfile = '../data/2018_10/EarlyV2_GM_LIN_unforced_damped_restart_moorings.mat';

mooringsfile = '/Volumes/Samsung_T5/nsf_iwv/2018_11/EarlyV2_GMexp_NL_forced_damped_64cube_moorings.mat';

load(mooringsfile)

nT = length(t);
nDepths = length(depths);

dt = t(2)-t(1);
S = zeros(nT+1,nDepths);
for iDepth = nDepths:-1:1
    cv_mooring = squeeze(u3d(:,iDepth,:) + sqrt(-1)*v3d(:,iDepth,:)).';
    [omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,[]);
    omega = [ -flipud(omega_p(2:end)); omega_p];
    S(:,nDepths-iDepth+1) = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
end

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Horizontal Velocity Spectrum');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];


% GM = GarrettMunkSpectrumConstantStratification(N0/2,[-Lz 0],latitude);
S_theory = GM.HorizontalVelocitySpectrumAtFrequencies(flip(depths),omega)';
S_theory( S_theory<1e-4 ) = nan;
semilogy(omega*86400/(2*pi),S_theory,'--','LineWidth', 2)

hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega*86400/(2*pi),S, 'LineWidth', 2)

hold on, 

xlabel('frequency (cycles per day)'), ylabel('m^2/s^2'), title('horizontal velocity power spectrum')

labels = cell(length(depths),1);
for iDepth=length(depths):-1:1
   labels{iDepth} = sprintf('%d m',round(depths(nDepths-iDepth+1))); 
end
% for i=(length(depths)+1):length(depths)*2
%    labels{i} = sprintf('%d m, GM',round(depths(i-3))); 
% end
legend(labels)

xlim([-12 12])
ylim([1e-3 2e2])

print('-depsc','HorizontalVelocitySpectrum.eps')