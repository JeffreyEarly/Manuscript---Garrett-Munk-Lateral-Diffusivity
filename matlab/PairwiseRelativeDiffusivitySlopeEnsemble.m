scaleFactor = 1;
LoadFigureDefaults;

runtype = 'nonlinear';

if strcmp(runtype,'linear')
    load('../data/2019_01/particles_LIN.mat');
    thetitle = 'Pairwise particle dispersion, linear simulation';
elseif strcmp(runtype,'nonlinear')
    load('../data/2019_01/particles_NL.mat');
    thetitle = 'Pairwise particle dispersion, nonlinear simulation';
else
    error('invalid run type.');
end

nLevels = size(x,2)/floatsPerLevel;

zLevel = 5;
zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
edges = CreateBinEdgesForInitialSeparation(t,x(:,zLevelIndices),y(:,zLevelIndices));

[r0,D2] = PairwiseMeanSquareSeparation( t, x(:,zLevelIndices), y(:,zLevelIndices), edges );
kappa_r = zeros(size(r0));
kappa_r_err = zeros(size(r0));
for iBin = 1:size(kappa_r,2)
    [slope, slope_err] = LinearLeastSquaresFit(t,D2(:,iBin),1);
    kappa_r(iBin) = slope/4;
    kappa_r_err(iBin) = slope_err/4;
end

kappa_r(isnan(r0)) = [];
kappa_r_err(isnan(r0)) = [];
r0(isnan(r0)) = [];

% The pairwise diffusivity is double the usual diffusivity
kappa = kappa_r/2;
kappa_err = kappa_r_err/2;

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','PairwiseParticleDisperion');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

subplot(2,1,1)
plot(t/86400,D2/1e6)
xlabel('days', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km^2', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title(thetitle, 'FontSize', figure_axis_label_size, 'FontName', figure_font);
set( gca, 'FontSize', figure_axis_tick_size);

subplot(2,1,2)
errorbar(r0/1e3, kappa,2*kappa_err)
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('m^2/s', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylim([0 1.1*max(kappa+2*kappa_err)])
set( gca, 'FontSize', figure_axis_tick_size);

print('-dpng', '-r300', sprintf('../figures/PairwiseParticleDispersion-%s.png',runtype))