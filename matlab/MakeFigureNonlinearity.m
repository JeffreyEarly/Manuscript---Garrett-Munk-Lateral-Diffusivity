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
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2019_05/';
end

if strcmp(runtype,'linear')

elseif strcmp(runtype,'nonlinear')
     load('/Volumes/Samsung_T5/nsf_iwv/2019_05//EarlyV2_GM_NL_forced_damped_restart_decomp.mat');
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end

% we need this to compute the line of constant frequency
WM = WintersModel(file);
wavemodel = WM.wavemodel;


% These are the values iMode, iK used in the MakeFiguresAutocorrelations
ExampleWaveJ = [8 5];
ExampleWaveK = [8 55];

ExampleVortexJ = 8;
ExampleVortexK = 55;

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
R = exp(2*(lambda_x+lambda_z)*(t(end)-t(1)));

% Let's contour the area that retains 90% of its value
C = contourc(j_diss,k_diss,R,0.90*[1 1]);
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
% Wave nonlinearity
%

FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

jpcolor( k, j, (waveHKEFromVariance./waveHKE)'); xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

scatter( k(ExampleWaveK), ExampleWaveJ, 25, 0*[1 1 1], 'filled' );
scatter( k(ExampleWaveK), ExampleWaveJ, 8, 1*[1 1 1], 'filled' );

[~,~,~,kk] = wavemodel.internalModes.ModesAtFrequency(2*2*pi/86400);
plot( kk,1:length(kk), '--', 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])

title('Wave Nonlinearity', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

yticks(ticks_y)
yticklabels(labels_y)
ylabel('vertical mode number', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

colormap( cmocean('dense') );
cb = colorbar('eastoutside');
caxis([0 1])

print('-dpng', '-r300', sprintf('../figures_2019_05/WaveNonlinearity.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vortex nonlinearity
%

FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

jpcolor( k, j, (vortexHKEFromVariance./vortexHKE)'); xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

scatter( k(ExampleVortexK), ExampleVortexJ, 25, 0*[1 1 1], 'filled' );
scatter( k(ExampleVortexK), ExampleVortexJ, 8, 1*[1 1 1], 'filled' );

title('Vortex Nonlinearity', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

yticks(ticks_y)
yticklabels(labels_y)
ylabel('vertical mode number', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

colormap( cmocean('dense') );
cb = colorbar('eastoutside');
caxis([0 1])

print('-dpng', '-r300', sprintf('../figures_2019_05/VortexNonlinearity.png'))