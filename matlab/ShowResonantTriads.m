dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown';

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

return

shouldShift = 1;
if shouldShift == 1
    shift = @(a) fftshift(a);
else
    shift = @(a) a;
end

k = shift(wavemodel.k);
l = shift(wavemodel.l);

K = shift(wavemodel.K(:,:,1));
L = shift(wavemodel.L(:,:,1));
Kh = shift(wavemodel.Kh(:,:,1));
Omega = shift(wavemodel.Omega(:,:,1)/wavemodel.f0);

RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);
RedundantCoefficients2D = shift(RedundantCoefficients(:,:,1));

% individual wave I want to find triads for...
iK = 70;
iL = 17;

k_i = K(iK,iL);
l_i = L(iK,iL);
omega_i = Omega(iK,iL);

% These are all possible combinations (including self and redundants) of
% resultant vectors.
k_p = k_i + K;
l_p = l_i + L;
omega_p = omega_i + Omega;

% let's make these column vectors
k_p = reshape(k_p,[],1);
l_p = reshape(l_p,[],1);
omega_p = reshape(omega_p,[],1);

% first, find the location of the resultant vector (might be nan!)
k_bin = discretize(k_p,k);
l_bin = discretize(l_p,l);

% find the resultant vectors that are *not* nan...
omega_bin = nan(size(k_bin));
validIndices = find(~isnan(k_bin) & ~isnan(l_bin));
linearIndices = sub2ind(size(Omega),k_bin(validIndices),l_bin(validIndices));
%...and find the associated frequencies
omega_bin(validIndices) = Omega(linearIndices);

% 

return;

Kh = Kh(~RedundantCoefficients2D);
Omega = Omega(~RedundantCoefficients2D);

[k,ia,ic] = unique(Kh);
omega = Omega(ia);

[k,I] = sort(k);
omega = abs(omega(I));

% this will be our tolerance... maybe median(diff(k)) is better? It's much
% smaller
dk = k(2)-k(1);
domega = omega(2)-omega(1);

% individual wave I want to find triads for...
iK = 10;
ki = k(iK);
omegai = omega(iK);

% add k, add omega
kp = abs(ki - k(iK:end));
omegap = abs(omegai - omega(iK:end));

% discreteize find the nearest bin for each point, given some bin edges. So
% we try to find the nearest k and nearest omega. If things are resonant,
% the two bins should be close to each other (or the same!).
k_bin = discretize(kp,k);
omega_bin = discretize(omegap,omega);

figure, plot(abs(k_bin-omega_bin))

% resonant_condition = abs(k_bin-omega_bin)<2;

