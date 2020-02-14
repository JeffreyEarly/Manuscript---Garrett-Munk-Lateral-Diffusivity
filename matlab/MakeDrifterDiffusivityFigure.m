scaleFactor = 1;
LoadFigureDefaults;

if 1
    load('/Users/jearly/Documents/ProjectRepositories/LatMix/drifters/observations/griddedRho1DriftersUnstrained.mat');
    drifterIndices = 1:9;
    
    x = x(:,drifterIndices);
    y = y(:,drifterIndices);
    titleText = 'Latmix Site 1 drifter dispersion, 6 days';
    outputFile = 'DispersionLatmix.eps';
else
    % load('/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_restart_particles.mat');
    load('/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM10_128x128x129_particles.mat');
    
    % spacing is 21
    drifterIndices = [3 24 43 44 45 46 47 66 87];
    
    x = x(1:1057,drifterIndices);
    y = y(1:1057,drifterIndices);
    titleText = 'Internal wave model float dispersion, 6 days';
    outputFile = 'DispersionModel.eps';
end





[x_com, y_com, q, r] = CenterOfMass( x, y );

x_sub = q;
y_sub = r;

q = q - q(1,:);
r = r - r(1,:);

figure, plot(q,r)
set( gca, 'FontSize', figure_axis_tick_size);
title(titleText, 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xlabel('m', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('m', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
axis equal

print('-depsc',sprintf('../figures/%s',outputFile))