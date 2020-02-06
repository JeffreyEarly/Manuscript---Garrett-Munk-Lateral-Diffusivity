GMEnergyLevel = 10^(0.5);
fileWithOutExtension = sprintf('/Users/jearly/Documents/nsf_iwv/DiffusivityExperiment_GM%02d_256x256x129',round(GMEnergyLevel*10));
file = sprintf('%s.nc',fileWithOutExtension);

t=ncread(file,'t');
x=ncread(file,'x-position').';
y=ncread(file,'y-position').';
z=ncread(file,'z-position').';
floatsPerLevel = size(x,2)/ncreadatt(file, '/', 'nFloatLevels');

outputfile = sprintf('%s_particles.mat',fileWithOutExtension);
save(outputfile,'x','y','z','t','floatsPerLevel', 'file','GMEnergyLevel');
