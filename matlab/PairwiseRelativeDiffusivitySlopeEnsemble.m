load('../data/2018_12/particles_NL.mat');

tIndices = 1:30:length(t);

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
    [r2(zLevel,:),kappa_r(zLevel,:),std_error(zLevel,:)] = PairwiseRelativeDiffusivityVsDistance(t(tIndices),x_float,y_float,'powspec',theBins);
end

figure
for zLevel=1:nLevels
    xMean = sqrt(r2(zLevel,:));
    yMean = kappa_r(zLevel,:);
    yStdErr = std_error(zLevel,:);
    
    errorbar(xMean/1e3, yMean,2*yStdErr), hold on
    if (max(yMean+2*yStdErr)>maxDiffusivity)
        maxDiffusivity = max(yMean+2*yStdErr);
    end
end