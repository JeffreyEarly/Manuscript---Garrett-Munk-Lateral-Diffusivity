diffusivityMetric = 'relative';
diffusivityMethod = 'powspec';


maxDiffusivity = 1.0;

load('../data/2018_12/particles_LIN.mat');

tIndices = 1:4:length(t);

nLevels = size(x,2)/floatsPerLevel;

dx = abs(x(1,2)-x(1,1));
theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
theBins(1) = 0;
nBins = length(theBins);

r2 = zeros(nLevels,nBins);
kappa_r = zeros(nLevels,nBins);
std_error = zeros(nLevels,nBins);
thelabels = cell(nLevels,1);

for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2(zLevel,:),kappa_r(zLevel,:),std_error(zLevel,:)] = PairwiseRelativeDiffusivityVsDistance(t(tIndices),x_float,y_float,diffusivityMethod,theBins);
end

figure
sp1 = subplot(2,1,1);
for zLevel=1:nLevels
    xMean = sqrt(r2(zLevel,:));
    yMean = kappa_r(zLevel,:);
    yStdErr = std_error(zLevel,:);
    
    yStdErr(xMean==0) = [];
    yMean(xMean==0) = [];
    xMean(xMean==0) = [];
    
    errorbar(xMean/1e3, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end

set(gca,'XTick',[])
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal *relative* diffusivity at different depths (%s method)',diffusivityMethod))
ylim([0 0.55])
text(250,0.9,'Linear')

load('../data/2018_12/particles_NL.mat');

tIndices = 1:4:length(t);

nLevels = size(x,2)/floatsPerLevel;

dx = abs(x(1,2)-x(1,1));
theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
theBins(1) = 0;
nBins = length(theBins);

r2 = zeros(nLevels,nBins);
kappa_r = zeros(nLevels,nBins);
std_error = zeros(nLevels,nBins);
thelabels = cell(nLevels,1);

for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2(zLevel,:),kappa_r(zLevel,:),std_error(zLevel,:)] = PairwiseRelativeDiffusivityVsDistance(t(tIndices),x_float,y_float,diffusivityMethod,theBins);
end

sp2 = subplot(2,1,2);
for zLevel=1:nLevels
    xMean = sqrt(r2(zLevel,:));
    yMean = kappa_r(zLevel,:);
    yStdErr = std_error(zLevel,:);
    
    yStdErr(xMean==0) = [];
    yMean(xMean==0) = [];
    xMean(xMean==0) = [];
    
    errorbar(xMean/1e3, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end

legend(thelabels,'Location', 'northwest')
text(250,0.9,'Nonlinear')

ylabel('diffusivity (m^2/s)')
xlabel('rms distance (km)')
ylim([0 1.0*maxDiffusivity])

subplot(sp1)
ylim([0 1.0*maxDiffusivity])

packfig(2,1)
print('-depsc2', sprintf('../figures/DiffusivityVsScaleLinNonlin_%s_%s.eps',diffusivityMetric,diffusivityMethod))

