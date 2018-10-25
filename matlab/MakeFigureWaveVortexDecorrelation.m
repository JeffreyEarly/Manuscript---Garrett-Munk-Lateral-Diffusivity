% Cim notes that the vortexHKE plot appears to show peaks along the axis of
% the radius of deformation for each mode.
%
% Check the individual autocorrelation functions leading to the structure
% in the vortex decorrelation time for the nonlinear run.
%
% In the region of wavenumber-mode space where the wave solutions do not
% decorrelated, how does the relative energy differ from a GM spectrum? In
% the decorrelated region, we don't expect this to match.

scaleFactor = 1;
LoadFigureDefaults;

runtype = 'linear';

% Needed just to pull out resolution information
if strcmp(runtype,'linear')

elseif strcmp(runtype,'nonlinear')
    load('../data/2018_10/EarlyV2_GM_NL_forced_damped_decomp.mat');
elseif strcmp(runtype,'nonlinear')
    load('../data/2018_10/EarlyV2_GM_NL_forced_damped_decomp.mat');
else
    error('invalid run type.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the damping scales
%

dx = x(2)-x(1);
dz = z(2)-z(1);
nu_x = (-1)^(p+1)*power(dx/pi,2*p) / T_diss;
nu_z = (-1)^(p+1)*power(dz/pi,2*p) / T_diss;

nK = length(k)/2 + 1;
k_diss = abs(k(1:nK));
j_diss = 0:max(j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(j)*dz))*J/2;

lambda_x = nu_x*(sqrt(-1)*K).^(2*p);
lambda_z = nu_z*(sqrt(-1)*M).^(2*p);
tau = max(t);
R = exp((lambda_x+lambda_z)*tau);

C = contourc(j_diss,k_diss,R,0.5*[1 1]);
n = C(2,1);
j_damp = C(1,1+1:n);
k_damp = C(2,1+1:n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the axis labels
%

ticks_x = [100;10];
labels_x = cell(length(ticks_x),1);
for i=1:length(ticks_x)
    labels_x{i} = sprintf('%d',round(ticks_x(i)));
end
ticks_x = 2*pi./(1e3*ticks_x);

ticks_y = [1;10;100];
labels_y = cell(length(ticks_y),1);
for i=1:length(ticks_y)
    labels_y{i} = sprintf('%d',round(ticks_y(i)));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decorrelation plots
%

colorAxis = reshape(waveDecorrelationTime,[],1);
colorAxis(isinf(colorAxis) | isnan(colorAxis)) = [];
colorAxisMin = min(colorAxis);
colorAxisMax = max(colorAxis);

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Wave-Vortex Decorrelation Time');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

p1 = subplot(1,2,1);
pcolor( k, j, (waveDecorrelationTime.')./86400 ), xlog, ylog, shading flat, hold on
% plot( sqrt(N2(GM.zInternal)), axis_depth, 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])
caxis([colorAxisMin colorAxisMax]/86400)
set( gca, 'FontSize', figure_axis_tick_size);
% set(gca, 'YTick', 1000*(-5:1:0));
% ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
title('wave decorrelation time', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

yticks(ticks_y)
yticklabels(labels_y)
ylabel('vertical mode number', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

p1.TickDir = 'out';

p2 = subplot(1,2,2);
pcolor( k, j, (vortexDecorrelationTime.')./86400 ), xlog, ylog, shading flat, hold on
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', []);
title('vortex decorrelation time', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

p2.TickDir = 'out';

colormap( cmocean('dense') );

cb = colorbar('eastoutside');
cb.Ticks = [0 2 4 6 8 10];
ylabel(cb, 'time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

print('-dpng', '-r300', sprintf('DecorrelationTime-%s.png',runtype))