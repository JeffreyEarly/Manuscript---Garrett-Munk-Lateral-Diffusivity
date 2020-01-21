% plot several examples of signal (u), with vertical axis set to std.
% deviation.

%fraction of HKE from variance to HKE

scaleFactor = 1;
LoadFigureDefaults;

runtype = 'nonlinear';
ReadOverNetwork = 0;
iFile = 2;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/';
end

energyLevel = [0.1; 1.0; 5.0];
if strcmp(runtype,'linear')
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_01xGM';
    files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart';
    files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_5xGM';
elseif strcmp(runtype,'nonlinear')
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM';
    files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart';
    files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM';
else
    error('invalid run type.');
end
figureFolder = '/Volumes/Samsung_T5/nsf_iwv/figures';

load(sprintf('%s_decomp.mat',files{iFile}));
file = files{iFile};


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

% Kraig Winters uses an e-fold time to set \nu_x in the hypervicous
% operator. We start by computing \nu_x and \nu_z.
dx = wavemodel.x(2)-wavemodel.x(1);
dz = wavemodel.z(2)-wavemodel.z(1);
nu_x = (-1)^(WM.p_x+1)*power(dx/pi,2*WM.p_x) / WM.T_diss_x;
nu_z = (-1)^(WM.p_z+1)*power(dz/pi,2*WM.p_z) / WM.T_diss_z;

nK = length(wavemodel.k)/2 + 1;
k_diss = abs(wavemodel.k(1:nK));
j_diss = 0:max(wavemodel.j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(wavemodel.j)*dz))*J/2;

lambda_x = nu_x*(sqrt(-1)*K).^(2*WM.p_x);
lambda_z = nu_z*(sqrt(-1)*M).^(2*WM.p_y);
tau = t(end)-t(1); tau = 86400/2;
R = exp(2*(lambda_x+lambda_z)*(tau));

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

% scatter( k(ExampleWaveK), ExampleWaveJ, 25, 0*[1 1 1], 'filled' );
% scatter( k(ExampleWaveK), ExampleWaveJ, 8, 1*[1 1 1], 'filled' );

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

print('-dpng', '-r300', sprintf('%s/WaveNonlinearity_%dxGM.png', figureFolder, round(energyLevel(iFile)*10) ))

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

% scatter( k(ExampleVortexK), ExampleVortexJ, 25, 0*[1 1 1], 'filled' );
% scatter( k(ExampleVortexK), ExampleVortexJ, 8, 1*[1 1 1], 'filled' );

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

print('-dpng', '-r300', sprintf('%s/VortexNonlinearity_%dxGM.png', figureFolder, round(energyLevel(iFile)*10) ))

% print('-dpng', '-r300', sprintf('figureFolder/VortexNonlinearity_%dxGM.png', round(energyLevel(iFile)*10) ))
% print('-dpng', '-r300', sprintf('../figures_2019_05/VortexNonlinearity.png'))