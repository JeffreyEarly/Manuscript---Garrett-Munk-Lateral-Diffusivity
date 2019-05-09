dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown';
decompfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown_decomp.nc';

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

k = ncread(decompfile, 'k');
l = ncread(decompfile, 'l');
j = ncread(decompfile, 'j');
t = ncread(decompfile, 't');


Ap = ncread(decompfile,'Ap_realp') + sqrt(-1)*ncread(decompfile,'Ap_imagp');
Am = ncread(decompfile,'Am_realp') + sqrt(-1)*ncread(decompfile,'Am_imagp');


[K,L,J] = ndgrid(k,l,j);
Kh = sqrt(K.*K + L.*L);
RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);

t_days = t/86400;
figure, plot( t_days, [squeeze(abs(Ap(9,2,1,:)).^2) squeeze(abs(Am(9,2,1,:)).^2) ])
return
decayrate = zeros(length(k),length(l));
iJ = 1;
for iK=1:length(k)
    for iL=1:length(l)
        energy = squeeze(abs(Ap(iK,iL,iJ,:)).^2)+squeeze(abs(Am(iK,iL,iJ,:)).^2);
        p = polyfit(t,energy,1);
        decayrate(iK,iL) = p(1);
    end
end


kn = (-length(k)/2):(length(k)/2 - 1);
ln = (-length(l)/2):(length(l)/2 - 1);

figure
jpcolor(kn,ln,fftshift(decayrate).')
colormap( cmocean('balance') );
caxis([-1e-10 1e-10])
cb = colorbar('eastoutside');



for iK=1:length(k)
    for iL=1:length(l)
        energy = squeeze(abs(Ap(iK,iL,iJ,:)).^2)+squeeze(abs(Am(iK,iL,iJ,:)).^2);
        decayrate(iK,iL) = (energy(end)-energy(1))/energy(1);
    end
end

figure
jpcolor(kn,ln,fftshift(decayrate).')
colormap( cmocean('balance') );
caxis([-1 1])
cb = colorbar('eastoutside');