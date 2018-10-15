% plot several examples of signal (u), with vertical axis set to std.
% deviation.

%fraction of HKE from variance to HKE

% scaleFactor = 1;
% LoadFigureDefaults;
% 
% 
% NonlinearSpindownFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
% NonlinearForcedFromInitialConditionsFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
% LinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';
% NonlinearSteadyStateFile = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';
% 
% 
% matfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyEtal_GM_NL_35e-11_36000s_restart_decomp.mat';
% load(matfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the damping scales
%

if ~exist('WM','var')
    WM = WintersModel(NonlinearSteadyStateFile);
    wavemodel = WM.wavemodel;
end

T_diss = 36000;
p = 3;
dx = wavemodel.x(2)-wavemodel.x(1);
dz = wavemodel.z(2)-wavemodel.z(1);
nu_x = (-1)^(p+1)*power(dx/pi,2*p) / T_diss;
nu_z = (-1)^(p+1)*power(dz/pi,2*p) / T_diss;

nK = length(wavemodel.k)/2 + 1;
k_diss = abs(wavemodel.k(1:nK));
j_diss = 0:max(wavemodel.j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(wavemodel.j)*dz))*J/2;

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
% Wave nonlinearity
%

FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

pcolor( k, j, (waveHKEFromVariance./waveHKE)'), xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

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

print('-dpng', '-r300', sprintf('WaveNonlinearity.png'))

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

pcolor( k, j, (vortexHKEFromVariance./vortexHKE)'), xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

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

print('-dpng', '-r300', sprintf('VortexNonlinearity.png'))