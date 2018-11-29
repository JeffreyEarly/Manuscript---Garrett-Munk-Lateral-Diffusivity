diffusivityMetric = 'relative';
diffusivityMethod = 'powspec';


maxDiffusivity = 1.0;

% load('../data/2018_11/particles_LIN.mat');

file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T154811_128x32x33.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T185207_128x32x33.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T193408_32x32x33.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T193717_64x64x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T195455_64x64x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T201915_64x32x33.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T202448_64x32x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T203604_64x32x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T204412_32x32x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-26T214847_64x64x65.nc';
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-27T065209_128x64x65.nc';

% nTidalWavelengths = 8
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-27T085340_512x128x65.nc';
% nTidalWavelengths = 4
file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-27T215253_512x128x65.nc';


t = ncread(file, 't');
nLevels = ncreadatt(file, '/', 'nFloatLevels');
x = ncread(file, 'x-position')';
y = ncread(file, 'y-position')';
z = ncread(file, 'z-position')';
floatsPerLevel = size(x,2)/nLevels;

tIndices = 1:30:length(t);
tIndexMax = find(t>9*86400,1,'first');
% tIndexMax = length(t);
tIndices = 1:1:tIndexMax;

% tIndices = 1:length(t);

% nLevels = size(x,2)/floatsPerLevel;

dx = abs(x(1,2)-x(1,1));
theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
theBins(1) = 0;
nBins = length(theBins);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% distance = zeros(floatsPerLevel*(floatsPerLevel-1)/2,1);
% iteration = 1;
% for i=1:floatsPerLevel
%     for j=(i+1):floatsPerLevel
%         dx = x(1,i) - x(1,j);
%         dy = y(1,i) - y(1,j);
%         distance(iteration) = sqrt( dx*dx + dy*dy );
%         iteration = iteration+1;
%     end
% end
% 
% dx = max(distance)/sqrt(floatsPerLevel);
% theBins = dx*(1:sqrt(floatsPerLevel)) + dx/2;
% theBins(1) = 0;
% nBins = length(theBins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r2 = zeros(nLevels,nBins);
kappa_r = zeros(nLevels,nBins);
std_error = zeros(nLevels,nBins);
thelabels = cell(nLevels,1);

for zLevel=1:nLevels
    zLevelIndices = (zLevel-1)*floatsPerLevel + (1:10:floatsPerLevel);
    
    x_float = x(tIndices,zLevelIndices);
    y_float = y(tIndices,zLevelIndices);
    z_float = z(tIndices,zLevelIndices);
    
    thelabels{zLevel} = sprintf('%d meters',round(mean(z_float(1,:))));
    [r2(zLevel,:),kappa_r(zLevel,:),std_error(zLevel,:)] = PairwiseRelativeDiffusivityVsDistance(t(tIndices),x_float,y_float,diffusivityMethod,theBins);
end

% t = t(tIndices);
% for i=1:2
%     for j=(i+1):2
%         du = diff(x_float(:,i))./diff(t) - diff(x_float(:,j))./diff(t);
%     end
% end

figure
% sp1 = subplot(2,1,1);
for zLevel=1:nLevels
    xMean = sqrt(r2(zLevel,:));
    yMean = kappa_r(zLevel,:);
    yStdErr = std_error(zLevel,:);
    
    errorbar(xMean/1e3, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end

xlabel('rms distance (km)')
ylabel('diffusivity (m^2/s)')
title(sprintf('horizontal *relative* diffusivity at different depths (%s method)',diffusivityMethod))
ylim([0 1.1*maxDiffusivity])

return

load('../data/2018_11/particles_NL.mat');


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
print('-depsc2', sprintf('DiffusivityVsScaleLinNonlin_%s_%s.eps',diffusivityMetric,diffusivityMethod))

