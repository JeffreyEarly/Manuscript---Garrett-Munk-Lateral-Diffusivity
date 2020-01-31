runtype = 'nonlinear';
ReadOverNetwork = 0;

figure('Name','LateralDiffusivityOfParticles-vs-Energy')
for i=1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Read in the nonlinear files, determine their actual energies, then go 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == 1
        files{1} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_01xGM';
        files{2} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_03xGM';
        files{3} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_restart';
        files{4} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_2xGM';
        files{5} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_5xGM';
        nFiles = 5;
        
        energyLevelMeasured = zeros(nFiles,1);
        for iFile=1:nFiles
            WM = WintersModel(files{iFile});
            wavemodel = WM.wavemodel;
            [t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
            wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
            energyLevelMeasured(iFile) = sum( abs(wavemodel.Amp_minus(:)).^2 + abs(wavemodel.Amp_plus(:)).^2 )/E;
        end
        
        energyLevel = energyLevelMeasured;
        files{1} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_01xGM_particles.mat';
        files{2} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_03xGM_particles.mat';
        files{3} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_restart_particles.mat';
        files{4} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_2xGM_particles.mat';
        files{5} = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_5xGM_particles.mat';
        nFiles = 5;
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read in the linear files
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        runtype = 'linear';
        energyLevel = 10.^[-1.0; -0.5; 0; 0.5];
        files{1} = '/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM01_128x128x129_particles.mat';
        files{2} = '/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM03_128x128x129_particles.mat';
        files{3} = '/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM10_128x128x129_particles.mat';
        files{4} = '/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM32_128x128x129_particles.mat';
        files{5} = '/Volumes/Samsung_T5/nsf_iwv/EarlyLinear/DiffusivityExperiment_GM100_128x128x129_particles.mat';
        nFiles = 4;
    end
    
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
    m(i) = p(1)/mu(2);
    C(i) = exp(p(2)-p(1)*mu(1)/mu(2));
    
    
    errorbar(energyLevel,kappa_particles,3*kappa_err_particles,'LineWidth',2), hold on
    ax = gca;
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    plot(energyLevel,C(i)*energyLevel.^m(i),'--','LineWidth',2)
end
xlog,ylog
xlabel('energy level (GM)')
ylabel('\kappa (m^2/s)')
title(sprintf('Lateral Diffusivity at (%d km)^2',round(sqrt(mean(r2_particles))*1e-3)))
legend('nonlinear', sprintf('nonlinear fit kappa = %.2f GM^{%.2f}',C(1),m(1)),'linear', sprintf('linear fit kappa = %.2f GM^{%.2f}',C(2),m(2)),'Location','northwest') 
xlim([0.08 12])

% print('-depsc','LateralDiffusivityOfParticles-WintersNL-EarlyLin.eps')
