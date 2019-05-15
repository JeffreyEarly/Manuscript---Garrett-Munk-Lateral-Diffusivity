scaleFactor = 1;
LoadFigureDefaults;

mooringsfile = '../data/2019_05/EarlyV2_GM_NL_forced_damped_spunupHiTide_moorings.mat';
mooringsfile = '../data/2019_05/EarlyV2_GM_NL_forced_damped_moorings.mat';


% mooringsfile = '../data/2018_12/EarlyV2_GM_LIN_unforced_damped_restart_moorings.mat';

% mooringsfile = '/Volumes/Samsung_T5/nsf_iwv/2018_11/EarlyV2_GMexp_NL_forced_damped_64cube_moorings.mat';

load(mooringsfile)

nT = length(t);



nDepths = length(depths);
depthLoop = length(depths):-4:1;

dt = t(2)-t(1);
S = zeros(nT+1,length(depthLoop));
i = 1;
for iDepth = depthLoop
    cv_mooring = squeeze(u3d(:,iDepth,:) + sqrt(-1)*v3d(:,iDepth,:)).';
    [omega_p, Spp, Snn, Spn] = mspec(dt,cv_mooring,[]);
    omega = [ -flipud(omega_p(2:end)); omega_p];
    S(:,length(depthLoop)-i+1) = (1/(2*pi))*[flipud(vmean(Snn,2)); vmean(Spp(2:end,:),2)];
    i = i+1;
end

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Horizontal Velocity Spectrum');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];


GM = GarrettMunkSpectrumConstantStratification(N0,[-Lz 0],latitude);
S_theory = GM.HorizontalVelocitySpectrumAtFrequencies(flip(depths(depthLoop)),omega)';
S_theory( S_theory<1e-4 ) = nan;
semilogy(omega*86400/(2*pi),S_theory,'--','LineWidth', 2)

hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega*86400/(2*pi),S, 'LineWidth', 2)

hold on, 

xlabel('frequency (cycles per day)'), ylabel('m^2/s^2'), title('horizontal velocity power spectrum')

labels = cell(length(depthLoop),1);
i = 1;
for iDepth=depthLoop
   labels{i} = sprintf('%d m',round(depths(nDepths-iDepth+1))); 
   i = i+1;
end
% for i=(length(depths)+1):length(depths)*2
%    labels{i} = sprintf('%d m, GM',round(depths(i-3))); 
% end
legend(labels)

xlim([-12 12])
ylim([2e-2 4e2])

print('-depsc','../data/2019_05/figures_EarlyV2_GM_NL_forced_damped/HorizontalVelocitySpectrum.eps')

dOmega = omega(2)-omega(1);
E_uv_theory = sum(S_theory,1,'omitnan')*dOmega;
E_uv = sum(S,1)*dOmega;
GMRatio = (E_uv./E_uv_theory).'

figure, semilogy(omega*86400/(2*pi),S./S_theory)