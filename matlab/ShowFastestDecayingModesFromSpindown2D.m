scaleFactor = 1;
LoadFigureDefaults;

dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/2019_05/EarlyV2_GM_NL_forced_damped_spindown';
decompfile = '/Volumes/Samsung_T5/nsf_iwv/2019_05/EarlyV2_GM_NL_forced_damped_spindown_decomp.nc';

% dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindownPSI';
% decompfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindownPSI_decomp.nc';

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
RedundantCoefficients2D = RedundantCoefficients(:,:,1);
Omega = wavemodel.Omega/wavemodel.f0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure of total energy in a given mode (positive and negative wave
% energy)
%
% t_days = t/86400;
% figure, plot( t_days, [squeeze(abs(Ap(1,1,1,:)).^2) squeeze(abs(Am(1,1,1,:)).^2) ])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1D figure of absolute decay rate
%

% We're going to use the k axis (which is finer than the l axis) for these
% values.
Kh2D = Kh(:,:,1);

kAxis = k(1:length(k)/2);
kAxis = sort(unique(abs(Kh2D)));

dk = kAxis(2)-kAxis(1);
nJ = length(j);
omegaAxis = zeros(length(kAxis),nJ);
decayrate = zeros(length(kAxis),nJ);
nensemeble = zeros(length(kAxis),nJ);

for iJ=1:length(j)
    
    Omega2D = Omega(:,:,iJ);
    for iK=1:length(kAxis)
        indices = find( abs(Kh2D) > kAxis(iK)-dk/2 & abs(Kh2D) <= kAxis(iK)+dk/2 & RedundantCoefficients2D == 0 );
        energy = zeros(length(t),1);
        % add all the energies in this band together...
        for iIndex=1:length(indices)
            [ii,jj]=ind2sub(size(Kh2D),indices(iIndex));
            energy = energy + squeeze(abs(Ap(ii,jj,iJ,:)).^2)+squeeze(abs(Am(ii,jj,iJ,:)).^2);
            nensemeble(iK,iJ) = nensemeble(iK,iJ) + 1;
        end
        %...then compute slope.
        %     p = polyfit(t,energy,1);
        %     decayrate(iK) = p(1);
        decayrate(iK,iJ) = (energy(end)-energy(1))/(t(end)-t(1));
        omegaAxis(iK,iJ) = mean(abs(Omega2D(indices)));
    end
end

% decayrate_pcolor = zeros(size(decayrate)+1);
% decayrate_pcolor(2:end,2:end)=decayrate;
% kAxis_pcolor = cat(1,kAxis(1)-(kAxis(2)-kAxis(1)),kAxis);
% j_pcolor = cat(1,0,j);

kAxis_shift=kAxis(1:end-1)+diff(kAxis);
kAxis_shift(end+1) = kAxis(end);

kAxis_shift=kAxis+1e-1*kAxis(2);

figure
jpcolor(kAxis_shift,j,decayrate.'), shading flat
xlabel('wavenumber (radians/m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('mode', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
title('Wave mode decay rate (100 days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xlog, ylog

n=25;
m = cmocean('balance');
m(end-n:end,:)=[];
m(1:n,:)=[];
colormap( m );
caxis([-1e-9 1e-9])
cb = colorbar('eastoutside');
cb.Label.String = 'energy rate (m^3/s^3)';

hold on
for fn=([1.15 2 3 4])
    [~,~,~,kk] = wavemodel.internalModes.ModesAtFrequency(fn*2*pi/86400);
    scatter( kk,(1:length(kk))+0.5, 36, 0*[1 1 1], 'filled' );
    scatter( kk,(1:length(kk))+0.5, 25, 1*[1 1 1], 'filled' );
end
% scatter( kk,1:length(kk), '--', 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])

% plot( kk,1:length(kk), '--', 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])



print('-dpng', '-r300', '../figures/energy_loss_during_spindown2D.png')

return

totalenergy = squeeze(sum(sum(abs(Ap(:,:,iJ,:)).^2,1),2))+squeeze(sum(sum(abs(Am(:,:,iJ,:)).^2,1),2));
p = polyfit(t,totalenergy,1);

negIndex = find( decayrate < 0,1,'first');
fprintf('For j=%d, the first negative decay rate occurs at %.2f f_0\n',iJ,omegaAxis(negIndex));

xlimits = [1 max(omegaAxis)];
xlimits = [1 10];
ylimits = [-1.6e-9 1.6e-9];

figure('Name',sprintf('DeltaEnergy, j=%d',iJ))
subplot(2,1,1)
% plot(omegaAxis,p(1)*ones(size(omegaAxis)),'Color',0.1*[1 1 1]), hold on
plot(omegaAxis,zeros(size(omegaAxis)),'Color',0.5*[1 1 1]), hold on
ax = gca;
ax.ColorOrderIndex = 1;
ax.XGrid = 'on';
plot(omegaAxis,decayrate,'LineWidth',2)
vlines(omegaAxis(negIndex),'g--')
ylabel('energy rate')
xlim(xlimits)
if iJ==1
   ylim(ylimits) 
end
subplot(2,1,2)
plot(omegaAxis,nensemeble,'LineWidth',2)
ylabel('number of modes')
xlabel('frequency (f_0)')
xlim(xlimits)
if iJ==1
   ylim([0 35]) 
end
packfig(2,1)

% print('-depsc2', '../figures/energy_loss_during_spindownPSI_zoomed.eps')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Absolute decay rate
%

decayrate = zeros(length(k),length(l));
iJ = 1;
for iK=1:length(k)
    for iL=1:length(l)
        energy = squeeze(abs(Ap(iK,iL,iJ,:)).^2)+squeeze(abs(Am(iK,iL,iJ,:)).^2);
        p = polyfit(t,energy,1);
        decayrate(iK,iL) = p(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2D figure of absolute decay rate
%

kn = (-length(k)/2):(length(k)/2 - 1);
ln = (-length(l)/2):(length(l)/2 - 1);

figure
jpcolor(kn,ln,fftshift(decayrate).')
colormap( cmocean('balance') );
caxis([-1e-10 1e-10])
cb = colorbar('eastoutside');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Relative decay rate
%

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