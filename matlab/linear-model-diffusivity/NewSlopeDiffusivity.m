file = '/Volumes/Samsung_T5/DiffusivityExperiment_2018-11-27T215253_512x128x65.nc';

t = ncread(file, 't');
nLevels = ncreadatt(file, '/', 'nFloatLevels');
x = ncread(file, 'x-position')';
y = ncread(file, 'y-position')';
z = ncread(file, 'z-position')';
floatsPerLevel = size(x,2)/nLevels;

zLevel = 1;
zLevelIndices = (zLevel-1)*floatsPerLevel + (1:10:floatsPerLevel);

x = x(:,zLevelIndices);
y = y(:,zLevelIndices);

nDrifters = size(x,2);

nPairs = nDrifters*(nDrifters-1)/2;
initialDistance = zeros(nPairs,1);
stride = 1;
iPair = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(1,iDrifter) - x(1,jDrifter);
        r = y(1,iDrifter) - y(1,jDrifter);
        
        initialDistance(iPair) = sqrt(q^2 + r^2);
        iPair = iPair + 1;
    end
end