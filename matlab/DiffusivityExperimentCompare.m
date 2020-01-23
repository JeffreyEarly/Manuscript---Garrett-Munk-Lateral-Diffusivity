wintersRunLinear = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart_particles.mat';
wintersRunNonlinear = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_particles.mat';
wavemodelRunSplineLowPass = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart_linrun_2020-01-21T184811_512x128x129.nc';
wavemodelRunSplineNoVortex = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart_linrun_2020-01-21T185216_512x128x129.nc';
nFiles = 4;

wintersRunLinear = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_01xGM_particles.mat';
wintersRunNonlinear = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_particles.mat';
wavemodelRunSplineLowPass = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_linrun_2020-01-22T114847_512x128x129.nc';
nFiles = 3;

file = wavemodelRunSplineLowPass;
spline_lowpass.t=ncread(file,'t');
spline_lowpass.x=ncread(file,'x-position').';
spline_lowpass.y=ncread(file,'y-position').';
spline_lowpass.z=ncread(file,'z-position').';

file = wavemodelRunSplineNoVortex;
spline_novortex.t=ncread(file,'t');
spline_novortex.x=ncread(file,'x-position').';
spline_novortex.y=ncread(file,'y-position').';
spline_novortex.z=ncread(file,'z-position').';

% load(wintersRun);
% 
% t1 = t_wm_lin;
% x1 = x_wm_lin;
% y1 = y_wm_lin;
% z1 = z_wm_lin;
% 
% % t1 = t;
% % x1 = x;
% % y1 = y;
% % z1 = z;
% 
% t2 = t_wm_spl;
% x2 = x_wm_spl;
% y2 = y_wm_spl;
% z2 = z_wm_spl;
% % 
% % t2 = t;
% % x2 = x;
% % y2 = y;
% % z2 = z;
% 
% validIndices = 1:(min(length(t1),length(t2))-1);
% 
% dx = x1(validIndices,:)-x2(validIndices,:);
% dy = y1(validIndices,:)-y2(validIndices,:);
% dz = z1(validIndices,:)-z2(validIndices,:);
% 
% figure
% plot(x1(validIndices,1:10),y1(validIndices,1:10)), hold on
% plot(x2(validIndices,1:10),y2(validIndices,1:10))


t_particles = cell(nFiles,1);
D2_particles = cell(nFiles,1);
r2_particles = zeros(nFiles,1);
for iFile=1:nFiles
    if iFile == 1
        load(wintersRunLinear);
    elseif iFile == 2
        load(wintersRunNonlinear);
    elseif iFile == 3
        t = spline_lowpass.t(1:length(spline_lowpass.t)-1,:);
        x = spline_lowpass.x(1:length(spline_lowpass.t)-1,:);
        y = spline_lowpass.y(1:length(spline_lowpass.t)-1,:);
    elseif iFile == 4
        t = spline_novortex.t(1:length(spline_novortex.t)-1,:);
        x = spline_novortex.x(1:length(spline_novortex.t)-1,:);
        y = spline_novortex.y(1:length(spline_novortex.t)-1,:);
    end
   
   t_particles{iFile} = t;
   D2_particles{iFile} = zeros(size(t));
   
   nLevels = 5;
   for zLevel = 1:nLevels
       nLevels = size(x,2)/floatsPerLevel;
       zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
       
       x_float = x(:,zLevelIndices);
       y_float = y(:,zLevelIndices);
       
       % The bin edge of 30km is chosen so that the rms separation of the
       % particles matches the dye.
       [D2_level,r2] = PairwiseRelativeDispersion( t, x_float, y_float, [0 30e3 Inf] );
       
       D2_particles{iFile} = D2_particles{iFile} + D2_level(:,1);
   end
   
   r2_particles(iFile) = r2(1);
   
   % one factor of 2 to average m_xx and m_yy, another factor to convert
   % from relative diffusivity
   D2_particles{iFile} = (D2_particles{iFile}/nLevels)/4;
end

figure('Name','LateralDispersion-vs-Energy')
for iFile=1:nFiles
   plot(t_particles{iFile}/86400,D2_particles{iFile}*1e-6,'LineWidth',2), hold on
end
title(sprintf('Lateral Dispersion at (%d km)^2',round(sqrt(mean(r2_particles))*1e-3)))
legend('winters linear', 'winters nonlinear','linear model low pass') %, 'linear model no vortex'
ylabel('km^2')
xlabel('days')
% print('-depsc', 'LateralDispersion-vs-Lowpass.eps')

figure
iParticle = 1;
load(wintersRunLinear);
plot(x(1:length(spline_novortex.t)-1,iParticle),y(1:length(spline_novortex.t)-1,iParticle)), hold on
load(wintersRunNonlinear);
plot(x(1:length(spline_novortex.t)-1,iParticle),y(1:length(spline_novortex.t)-1,iParticle))
% plot(spline_novortex.x(1:length(spline_novortex.t)-1,iParticle),spline_novortex.y(1:length(spline_novortex.t)-1,iParticle))
plot(spline_lowpass.x(1:length(spline_lowpass.t)-1,iParticle),spline_lowpass.y(1:length(spline_lowpass.t)-1,iParticle))
legend('winters linear', 'winters nonlinear','linear model low pass') %, 'linear model no vortex'


