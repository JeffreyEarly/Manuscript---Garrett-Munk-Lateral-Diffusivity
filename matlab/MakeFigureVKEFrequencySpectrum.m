scaleFactor = 1;
LoadFigureDefaults;

mooringsfile = '../data/2018_10/EarlyV2_GM_NL_forced_damped_moorings.mat';
mooringsfile = '../data/2018_10/EarlyV2_GM_LIN_unforced_damped_moorings.mat';

load(mooringsfile)

nT = length(t);
nDepths = length(depths);

dt = t(2)-t(1);
S = zeros(round(nT/2),nDepths);
for iDepth = 1:nDepths
    v_mooring = squeeze(w3d(:,iDepth,:)).';
    [omega, Spp] = mspec(dt,v_mooring,[]);
    S(:,iDepth) = (1/(2*pi))*vmean(Spp,2);
end

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Vertical Velocity Spectrum');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];


GM = GarrettMunkSpectrumConstantStratification(N0,[-Lz 0],latitude);
S_theory = GM.VerticalVelocitySpectrumAtFrequencies(depths,omega)';
S_theory( S_theory<1e-4 ) = nan;
semilogy(omega*86400/(2*pi),S_theory,'--','LineWidth', 2)

hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega*86400/(2*pi),S, 'LineWidth', 2)

hold on, 

xlabel('frequency (cycles per day)'), ylabel('m^2/s^2'), title('vertical velocity power spectrum')

% labels = cell(length(depths)*2,1);
% for i=1:length(depths)
%    labels{i} = sprintf('%d m, model',round(depths(i))); 
% end
% for i=(length(depths)+1):length(depths)*2
%    labels{i} = sprintf('%d m, GM',round(depths(i-3))); 
% end
% legend(labels)

xlim([0 12])


% print('-depsc','HorizontalVelocitySpectrum2.eps')