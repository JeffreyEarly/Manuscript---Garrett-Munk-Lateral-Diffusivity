runtype = 'nonlinear';
ReadOverNetwork = 0;

energyLevel = [0.1; 1.0; 5.0];
files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_particles.mat';
files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_particles.mat';
files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM_particles.mat';

runtype = 'linear';
energyLevel = 10.^[-1.0; -0.5; 0; 0.5];
files{1} = '/Volumes/MoreStorage/DiffusivityExperiment_GM01_128x128x129_particles.mat';
files{2} = '/Volumes/MoreStorage/DiffusivityExperiment_GM03_128x128x129_particles.mat';
files{3} = '/Volumes/MoreStorage/DiffusivityExperiment_GM10_128x128x129_particles.mat';
files{4} = '/Volumes/MoreStorage/DiffusivityExperiment_GM32_128x128x129_particles.mat';
nFiles = 4;

t_particles = cell(nFiles,1);
D2_particles = cell(nFiles,1);
r2_particles = zeros(nFiles,1);
for iFile=1:nFiles
   load(files{iFile});
   
   t_particles{iFile} = t;
   D2_particles{iFile} = zeros(size(t));
   
   nLevels = 5;
   for zLevel = 1:nLevels
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

kappa_particles = zeros(size(energyLevel));
kappa_err_particles = zeros(nFiles,1);
for iFile=1:nFiles
   % factor of two from the defition
   [D2_coeff,D2_err] = LinearLeastSquaresFit(t_particles{iFile},D2_particles{iFile});
   kappa_particles(iFile) = D2_coeff(2)/2;
   kappa_err_particles(iFile) = D2_err(2)/2;
end


[p,S,mu]=polyfit(log(energyLevel),log(kappa_particles),1);
m = p(1)/mu(2);
C = exp(p(2)-p(1)*mu(1)/mu(2));

figure('Name','LateralDiffusivityOfParticles-vs-Energy')
plot(energyLevel,C*energyLevel.^m), hold on
scatter(energyLevel,kappa_particles),xlog,ylog
xlabel('energy level (GM)')
ylabel('\kappa (m^2/s)')
title(sprintf('kappa = %.2f GM^{%.2f}',C,m));