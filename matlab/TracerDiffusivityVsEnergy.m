runtype = 'nonlinear';
ReadOverNetwork = 0;

energyLevel = [0.1; 1.0; 5.0];
files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_01xGM_tracer_patch.mat';
files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_restart_tracer_patch.mat';
files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_5xGM_tracer_patch.mat';

figure('Name','LateralDiffusivityOfTracer')

kappa = zeros(size(energyLevel));
for iFile=1:3
   load(files{iFile});
   D2 = (m_xx+m_yy)/2;
   [D2_coeff,D2_err] = LinearLeastSquaresFit(t,D2);
   kappa(iFile) = D2_coeff(2)/2;
   
   plot(t/86400,D2/1e6), hold on
   plot(t/86400,(D2_coeff(2)*t+D2_coeff(1))/1e6,'k--')
end

% log(kappa)=m*log(energyLevel)+b
% kappa=exp(b)*energyLevel^b
[p,S,mu]=polyfit(log(energyLevel),log(kappa),1);
m = p(1)/mu(2);
C = exp(p(2)-p(1)*mu(1)/mu(2));

figure('Name','LateralDiffusivityOTracer-vs-Energy')
plot(energyLevel,C*energyLevel.^m), hold on
scatter(energyLevel,kappa),xlog,ylog
xlabel('energy level (GM)')
ylabel('\kappa (m^2/s)')
title(sprintf('kappa = %.2f GM^{%.2f}',C,m));