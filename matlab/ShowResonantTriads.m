% Resonant Triads
%
% The approach taken here is to find the *most* resonant triad for each
% given wavenumber. For a forced simple harmonic oscillator, the amplitude
% goes like 1/(omega-omega_0) where omega_0 is the natural frequency. So we
% use value to determine 'most'..

dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown';

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

% Discretize doesn't work unless you shift.
shouldShift = 1;
if shouldShift == 1
    shift = @(a) fftshift(a);
else
    shift = @(a) a;
end

nModes = 3;
modes = 1:3;

k = shift(wavemodel.k);
l = shift(wavemodel.l);

dk = k(2)-k(1);
dl = l(2)-l(1);

K = shift(wavemodel.K(:,:,1));
L = shift(wavemodel.L(:,:,1));
M = shift(wavemodel.M(:,:,1));
Kh = shift(wavemodel.Kh(:,:,1));
Omega = shift(wavemodel.Omega(:,:,1)/wavemodel.f0);

RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);
RedundantCoefficients2D = shift(RedundantCoefficients(:,:,1));

% Create a structure to log the most resonant triad for each wavenumber
nothing = NaN(size(K));
resonance = struct('i2',nothing,'j2',nothing,'k2',nothing,'l2',nothing,'omega2',nothing,'i3',nothing,'j3',nothing,'k3',nothing,'l3',nothing,'omega3',nothing,'Amp',nothing);

% individual wave I want to find triads for...
iK = 70;
iL = 17;

for iK=1:length(k)
    for iL=1:length(l)        
%         if RedundantCoefficients2D(iK,iL)==1
%             continue;
%         end
        
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
        
        Amp = 1./abs(omega_bin-omega_p);
        Amp(logical(reshape(RedundantCoefficients2D,[],1))) = NaN;
        [maxAmp,maxAmpIndex] = max(Amp,[],'omitnan');
        
        [i,j] = ind2sub(size(K),maxAmpIndex);
        resonance.i2(iK,iL) = i;
        resonance.j2(iK,iL) = j;
        resonance.k2(iK,iL) = K(i,j);
        resonance.l2(iK,iL) = L(i,j);
        resonance.omega2(iK,iL) = Omega(i,j);
        
        resonance.i3(iK,iL) = k_bin(maxAmpIndex);
        resonance.j3(iK,iL) = l_bin(maxAmpIndex);
        resonance.k3(iK,iL) = K(k_bin(maxAmpIndex),l_bin(maxAmpIndex));
        resonance.l3(iK,iL) = L(k_bin(maxAmpIndex),l_bin(maxAmpIndex));
        resonance.omega3(iK,iL) = omega_bin(maxAmpIndex);
        
        resonance.omega12(iK,iL) = omega_p(maxAmpIndex);
        
        resonance.Amp(iK,iL) = maxAmp;
    end
end

[sortedResonance,sortOrder] = sort(reshape(resonance.Amp,[],1),'descend');
originalIndex = (1:length(sortedResonance))';
originalIndex = originalIndex(sortOrder);

for iResonant=1:3
    maxAmp = sortedResonance(iResonant);
    maxAmpIndex = originalIndex(iResonant);
    
    [i,j] = ind2sub(size(K),maxAmpIndex);
    fprintf('\nResonant pair %d is (i,j)=(%d,%d) added to (i,j)=(%d,%d) to get (i,j)=(%d,%d)\n',iResonant,i,j,resonance.i2(i,j),resonance.j2(i,j),resonance.i3(i,j),resonance.j3(i,j));
    fprintf('This is (k,l,omega): (%.2g, %.2g,%.2f f0) + (%.2g, %.2g,%.2f f0) = (%.2g, %.2g,%.2f f0)\n',K(i,j),L(i,j),Omega(i,j),resonance.k2(i,j),resonance.l2(i,j),resonance.omega2(i,j),resonance.k3(i,j),resonance.l3(i,j),resonance.omega3(i,j));
end

figure
jpcolor(k,l,resonance.Amp');
% hold on
% scatter(k_i+dk/2,l_i+dl/2,25,0*[1 1 1],'filled')

return



figure
jpcolor(k,l,RedundantCoefficients2D.')
hold on
scatter(k_i+dk/2,l_i+dl/2,25,0*[1 1 1],'filled')

figure
jpcolor(k,l,Omega.')
hold on
scatter(k_i+dk/2,l_i+dl/2,25,0*[1 1 1],'filled')


% figure
% plot(1./abs(omega_bin-omega_p))

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

