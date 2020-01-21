runtype = 'nonlinear';
ReadOverNetwork = 0;

energyLevel = [0.1; 1.0; 5.0];
files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_particles.mat';
files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_particles.mat';
files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM_particles.mat';



kappa = zeros(size(energyLevel));
for iFile=1:3
   load(files{iFile});
   
   zLevel = 3;
   nLevels = size(x,2)/floatsPerLevel;
   zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
   tIndices = 1:4:length(t);
   
   t_float = t(tIndices);
   x_float = x(tIndices,zLevelIndices);
   y_float = y(tIndices,zLevelIndices);
   
%    D2 = mean((x_float-x_float(1,:)).^2 + (y_float-y_float(1,:)).^2,2);
%    kappa(iFile) = LinearLeastSquaresFit(t_float,D2,1)/4;
   
%    D2 = mean((x-x(1,:)).^2 + (y-y(1,:)).^2,2);
%    kappa(iFile) = LinearLeastSquaresFit(t,D2,1)/4;
   
%    [r2_all, kappa_r_all] = PairwiseRelativeDiffusivity(t_float, x_float, y_float, 'powspec');
%   kappa(iFile) = mean(kappa_r_all)/2;
   [r2_all, kappa_r_all] = PairwiseRelativeDiffusivityFromSlope(t, x, y, [0 1685] );
   kappa(iFile) = kappa_r_all(1)/2;
   
end

% log(kappa)=m*log(energyLevel)+b
% kappa=exp(b)*energyLevel^b
[p,S,mu]=polyfit(log(energyLevel),log(kappa),1);
m = p(1)/mu(2);
C = exp(p(2)-p(1)*mu(1)/mu(2));

figure('Name','LateralDiffusivityOfParticles-vs-Energy')
plot(energyLevel,C*energyLevel.^m), hold on
scatter(energyLevel,kappa),xlog,ylog
xlabel('energy level (GM)')
ylabel('\kappa (m^2/s)')
title(sprintf('kappa = %.2f GM^{%.2f}',C,m));