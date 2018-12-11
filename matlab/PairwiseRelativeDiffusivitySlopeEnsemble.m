load('../data/2018_12/particles_NL.mat');

nLevels = size(x,2)/floatsPerLevel;

zLevel = 5;
zLevelIndices = (zLevel-1)*floatsPerLevel + (1:floatsPerLevel);
    
x0 = x(1,zLevelIndices);
y0 = y(1,zLevelIndices);
edges = CreateBinEdgesForInitialSeparation(x0,y0);

[r0,D2] = PairwiseMeanSquareSeparation( t, x(:,zLevelIndices), y(:,zLevelIndices), edges );
kappa_r = zeros(size(r0));
kappa_r_err = zeros(size(r0));
for iBin = 1:size(kappa_r,2)
    [slope, slope_err] = LinearLeastSquaresFit(t,D2(:,iBin),1);
    kappa_r(iBin) = slope/4;
    kappa_r_err(iBin) = slope_err/4;
end

kappa_r(isnan(r0)) = [];
kappa_r_err(isnan(r0)) = [];
r0(isnan(r0)) = [];

figure,
subplot(2,1,1)
plot(t/86400,D2/1e6)
xlabel('days')
ylabel('km^2')
subplot(2,1,2)
errorbar(r0/1e3, kappa_r,2*kappa_r_err)
xlabel('km')
ylabel('m^2/s')
