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

runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end

WM = WintersModel(file);
wavemodel = WM.wavemodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Do the wave-vortex decomposition
%

maxFileIndex = WM.NumberOf3DOutputFiles;
t0 = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t');
t_max = WM.VariableFieldsFrom3DOutputFileAtIndex(maxFileIndex,'t');

% These first two lines actually do the decomposition
[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(maxFileIndex,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);

waveHKE_full = wavemodel.Ppm_HKE_factor .* ( abs(wavemodel.Amp_plus).^2 + abs(wavemodel.Amp_minus).^2 );
vortexHKE_full = wavemodel.P0_HKE_factor .* abs(wavemodel.B).^2;

% Now we need to figure out how to display this information...

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

% Thi is the final output axis for wavenumber
k = reshape(kAxis(1:(length(kAxis)-1)),[],1);

% Mode axis is just what we already have
j = wavemodel.j;

Kh = wavemodel.Kh;
RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);

nK = length(k);
nModes = length(j);

waveHKE = zeros(nK,nModes);
vortexHKE = zeros(nK,nModes);

for iMode = 1:1:nModes
    for iK = 1:1:nK
        indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(RedundantCoefficients(:,:,1)) );
        for iIndex = 1:length(indicesForK)
            [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
            waveHKE(iK,iMode) = waveHKE(iK,iMode) + waveHKE_full(i,m,iMode);
            vortexHKE(iK,iMode) = vortexHKE(iK,iMode) + vortexHKE_full(i,m,iMode);
        end
    end
end
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
R = exp(2*(lambda_x+lambda_z)*(t_max-t0));

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
% Wave/Vortex energy
%

% These are the values iMode, iK used in the MakeFiguresAutocorrelations
ExampleWaveJ = [8 5];
ExampleWaveK = [8 55];

ExampleVortexJ = 8;
ExampleVortexK = 55;

% Compute the Rossby radius of deformation for each geostrophic mode
% N0 = 5.2e-3;
% D = 4000;
% z = linspace(-D,0,1024)';
% im = InternalModesConstantStratification(N0,[-D 0],z,33);
% [F,G,h] = im.ModesAtFrequency(0.0);
% Lr = sqrt(9.81*h(j))/im.f0;
% kr = 1./Lr; % I think it should *not* be multiplied by 2*pi, because in the dispersion relation it adds directly to k and l.

FigureSize = [50 50 figure_width_2col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize,'Name','Wave-Vortex Energy');
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

p1 = subplot(1,3,1);
jpcolor( k, j, log10(waveHKE.') ); xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])
scatter( k(ExampleWaveK), ExampleWaveJ, 25, 0*[1 1 1], 'filled' );
scatter( k(ExampleWaveK), ExampleWaveJ, 8, 1*[1 1 1], 'filled' );
% plot( sqrt(N2(GM.zInternal)), axis_depth, 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])
caxis([-10 0])
set( gca, 'FontSize', figure_axis_tick_size);
% set(gca, 'YTick', 1000*(-5:1:0));
% ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
title('wave energy (m^3/s^2)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

yticks(ticks_y)
yticklabels(labels_y)
ylabel('vertical mode number', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

p1.TickDir = 'out';
cb1 = colorbar('eastoutside');
cb1.Ticks = [-10 -8 -6 -4 -2 0];
cb1.TickLabels = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}'};
% ylabel(cb, 'depth integrated energy (m^3/s^2)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

p2 = subplot(1,3,2);
jpcolor( k, j, log10(vortexHKE.') ); xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])
scatter( k(ExampleVortexK), ExampleVortexJ, 25, 0*[1 1 1], 'filled' );
scatter( k(ExampleVortexK), ExampleVortexJ, 8, 1*[1 1 1], 'filled' );
caxis([-10 -5])
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', []);
title('vortex energy (m^3/s^2)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

% hold on, plot(kr,j, 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])

[~,~,~,kk] = wavemodel.internalModes.ModesAtFrequency(2.0*2*pi/86400);
plot( kk,1:length(kk), '--', 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])

p2.TickDir = 'out';

colormap( cmocean('dense') );

cb2 = colorbar('eastoutside');
cb2.Ticks = [-10 -9 -8 -7 -6 -5];
cb2.TickLabels = {'10^{-10}', '10^{-9}', '10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}'};
% ylabel(cb2, 'm^3/s^2', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Energy fraction
%

HKE_fraction = waveHKE./(waveHKE+vortexHKE);

% FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
% fig1 = figure('Units', 'points', 'Position', FigureSize);
% set(gcf, 'Color', 'w');
% fig1.PaperUnits = 'points';
% fig1.PaperPosition = FigureSize;
% fig1.PaperSize = [FigureSize(3) FigureSize(4)];

p3 = subplot(1,3,3);
jpcolor( k, j, abs(log10(1-HKE_fraction.') )); xlog, ylog, shading flat, hold on
plot( k_damp, j_damp, 'LineWidth', 4, 'Color', 0*[1 1 1])
plot( k_damp, j_damp, 'LineWidth', 2, 'Color', [1 1 1])

% [~,~,~,kk] = im.ModesAtFrequency(2.0*2*pi/86400);
plot( kk,1:length(kk), '--', 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])

set( gca, 'FontSize', figure_axis_tick_size);

title('wave energy fraction', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

xticks(ticks_x)
xticklabels(labels_x)
xlabel('wavelength (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

% yticks(ticks_y)
% yticklabels(labels_y)
% ylabel('vertical mode number', 'FontSize', figure_axis_label_size, 'FontName', figure_font);


set(gca, 'YTick', []);
p3.TickDir = 'out';

colormap( cmocean('dense') );
cb3 = colorbar('eastoutside');
caxis(abs(log10(1-[(0.5) (0.99999)])))
cb3.Ticks = abs(log10(1-[(0.5) (0.9) (0.99) (0.999) (0.9999) (0.99999)]));
cb3.TickLabels = {'50%', '90%', '99%', '99.9%', '99.99%', '99.999%'};

p1.OuterPosition = [0 p1.OuterPosition(2) 0.28 p1.OuterPosition(4)];
p1.Position = [p1.Position(1) p1.Position(2) 0.20 p1.Position(4)];
p2.OuterPosition = [0.36 p2.OuterPosition(2) 0.28 p2.OuterPosition(4)];
p2.Position = [p2.Position(1) p2.Position(2) 0.20 p2.Position(4)];
p3.OuterPosition = [0.66 p3.OuterPosition(2) 0.28 p3.OuterPosition(4)];
p3.Position = [p3.Position(1) p3.Position(2) 0.20 p3.Position(4)];

print('-dpng', '-r300', sprintf('../figures/EnergyAndEnergyFraction-%s.png',runtype))
