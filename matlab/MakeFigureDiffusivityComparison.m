diffusivityMetric = 'relative';
diffusivityMethod = 'powspec';


maxDiffusivity = 0;

load('../data/2018_10/particles_LIN.mat');

tIndices = 1:length(t);

nLevels = size(x,2)/floatsPerLevel;

r2_r = cell(nLevels,1);
kappa_r = cell(nLevels,1);
kappa_r_corr = cell(nLevels,1);
r2_a = cell(nLevels,1);
kappa_a = cell(nLevels,1);
thelabels = cell(nLevels,1);
for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    if strcmp(diffusivityMetric,'relative')
        [r2_r{zLevel}, kappa_r{zLevel}, kappa_r_corr{zLevel}] = PairwiseRelativeDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
    elseif strcmp(diffusivityMetric,'absolute')
        [r2_a{zLevel}, kappa_a{zLevel}] = SingleParticleDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
    end
end

dx = abs(x(1,2)-x(1,1));
theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
theBins(1) = 0;


figure
sp1 = subplot(2,1,1);
for zLevel=1:nLevels
    if strcmp(diffusivityMetric,'relative')
        [xMean,yMean,yStdErr] = BinDataWithErrorBars(sqrt(r2_r{zLevel}),kappa_r{zLevel},theBins);
    elseif strcmp(diffusivityMetric,'absolute')
        [xMean,yMean,yStdErr] = BinDataWithErrorBars(sqrt(r2_a{zLevel}),kappa_a{zLevel},theBins);
    end
    errorbar(xMean, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end

set(gca,'XTick',[])
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal *relative* diffusivity at different depths (%s method)',diffusivityMethod))
ylim([0 0.55])

load('../data/2018_10/particles_NL.mat');

tIndices = 1:length(t);

n = floatsPerLevel;
nLevels = size(x,2)/floatsPerLevel;

r2_r = cell(nLevels,1);
kappa_r = cell(nLevels,1);
kappa_r_corr = cell(nLevels,1);
r2_a = cell(nLevels,1);
kappa_a = cell(nLevels,1);
theZlabels = cell(nLevels,1);
thelabels = cell(nLevels,1);
for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2_r{zLevel}, kappa_r{zLevel}, kappa_r_corr{zLevel}] = PairwiseRelativeDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
    [r2_a{zLevel}, kappa_a{zLevel}] = SingleParticleDiffusivity(t(tIndices),x_float,y_float,diffusivityMethod);
end

sp2 = subplot(2,1,2);
for zLevel=1:nLevels
    if strcmp(diffusivityMetric,'relative')
        [xMean,yMean,yStdErr] = BinDataWithErrorBars(sqrt(r2_r{zLevel}),kappa_r{zLevel},theBins);
    elseif strcmp(diffusivityMetric,'absolute')
        [xMean,yMean,yStdErr] = BinDataWithErrorBars(sqrt(r2_a{zLevel}),kappa_a{zLevel},theBins);
    end
    errorbar(xMean, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end

legend(thelabels,'Location', 'northwest')

ylabel('diffusivity (m^2/s)')
xlabel('distance (m)')
ylim([0 1.1*maxDiffusivity])

subplot(sp1)
ylim([0 1.1*maxDiffusivity])

packfig(2,1)
print('-depsc2', sprintf('DiffusivityVsScaleLinNonlin_%s_%s.eps',diffusivityMetric,diffusivityMethod))

