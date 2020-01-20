runtype = 'nonlinear';
ReadOverNetwork = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle dispersion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(runtype,'nonlinear') == 1
    energyLevel = [0.1; 1.0; 5.0];
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_particles.mat';
    files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_particles.mat';
    files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM_particles.mat';
    nFiles = 3;
else
    energyLevel = 1.0;
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart_particles.mat';
    nFiles = 1;
end

t_particles = cell(nFiles,1);
D2_particles = cell(nFiles,1);
r2_particles = zeros(nFiles,1);
kappa_particles = zeros(size(energyLevel));
for iFile=1:nFiles
   load(files{iFile});
   
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tracer dispersion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(runtype,'nonlinear') == 1
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_tracer_patch.mat';
    files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_tracer_patch.mat';
    files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM_tracer_patch.mat';
else
    energyLevel = 1.0;
    files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_LIN_unforced_damped_restart_tracer_patch.mat';
    nFiles = 1;
end

t_tracer = cell(nFiles,1);
D2_tracer = cell(nFiles,1);
r2_tracer = zeros(nFiles,1);
kappa_tracer = zeros(size(energyLevel));
for iFile=1:3
   load(files{iFile});
   
   t_tracer{iFile} = t;
   r2_tracer(iFile) = m_xx(1)+m_yy(1);
   
   % one factor of 2 to average m_xx and m_yy
   D2_tracer{iFile} = (m_xx+m_yy)/2 - (m_xx(1)+m_yy(1))/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Diffusivities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa_err_particles = zeros(nFiles,1);
kappa_err_tracer = zeros(nFiles,1);
for iFile=1:nFiles
   % factor of two from the defition
   [D2_coeff,D2_err] = LinearLeastSquaresFit(t_particles{iFile},D2_particles{iFile});
   kappa_particles(iFile) = D2_coeff(2)/2;
   kappa_err_particles(iFile) = D2_err(2)/2;
   
   [D2_coeff,D2_err] = LinearLeastSquaresFit(t_tracer{iFile},D2_tracer{iFile});
   kappa_tracer(iFile) = D2_coeff(2)/2;
end

% Relationship between diffusivity and energy
[p,S,mu]=polyfit([log(energyLevel); log(energyLevel)],[log(kappa_particles); log(kappa_tracer)],1);
m = p(1)/mu(2);
C = exp(p(2)-p(1)*mu(1)/mu(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dispersion Figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','LateralDispersion-vs-Energy')
for iFile=1:nFiles
   subplot(nFiles,1,iFile)
   plot(t_particles{nFiles-iFile+1}/86400,D2_particles{nFiles-iFile+1}*1e-6,'LineWidth',2), hold on
   plot(t_tracer{nFiles-iFile+1}/86400,D2_tracer{nFiles-iFile+1}*1e-6,'LineWidth',2)
   if iFile==1
      title(sprintf('Lateral Dispersion at (%d km)^2',round(sqrt(mean(r2_particles))*1e-3)))
      legend('particles','tracer','Location','northwest') 
   end
   ylabel('km^2')
end
xlabel('days')
packfig(nFiles,1)
print('-depsc', sprintf('../figures_2020_01/LateralDispersion-vs-Energy.eps'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Diffusivity Figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','LateralDiffusivity-vs-Energy')
errorbar(energyLevel,kappa_particles,3*kappa_err_particles,'LineWidth',2), hold on
errorbar(energyLevel,kappa_tracer,3*kappa_err_tracer,'LineWidth',2)
plot(energyLevel,C*energyLevel.^m,'k--','LineWidth',2)
xlog
ylog
xlabel('energy level (GM)')
ylabel('\kappa (m^2/s)')
title(sprintf('Lateral Diffusivity at (%d km)^2, kappa = %.2f GM^{%.2f}',round(sqrt(mean(r2_particles))*1e-3),C,m))
legend('particles','tracer','Location','northwest') 
print('-depsc', sprintf('../figures_2020_01/LateralDiffusivity-vs-Energy.eps'))

