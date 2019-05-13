% Resonant Triads
%
% The approach taken here is to find the *most* resonant triad for each
% given wavenumber. For a forced simple harmonic oscillator, the amplitude
% goes like 1/(omega-omega_0) where omega_0 is the natural frequency. So we
% use value to determine 'most'..

% If you walk through all coefficients (including negative k,l), then you
% want omega to always be positive.
%
% But then don't we need to consider the difference in omega?

dynamicalfile = '/Volumes/Samsung_T5/nsf_iwv/EarlyV2_GM_NL_forced_damped_spindown';

WM = WintersModel(dynamicalfile);
wavemodel = WM.wavemodel;

% 
normalizeByGMAmplitude = 0;


% Discretize doesn't work unless you shift.
shouldShift = 1;
if shouldShift == 1
    shift = @(a) fftshift(fftshift(a,1),2);
else
    shift = @(a) a;
end

nModes = 3;
modes = 1:nModes;

k = shift(wavemodel.k);
l = shift(wavemodel.l);
m = wavemodel.j(modes)*pi/wavemodel.Lz;

dk = k(2)-k(1);
dl = l(2)-l(1);
dm = m(2)-m(1);

K = shift(wavemodel.K(:,:,modes));
L = shift(wavemodel.L(:,:,modes));
M = shift(wavemodel.M(:,:,modes));
Kh = shift(wavemodel.Kh(:,:,modes));
Omega = abs(shift(wavemodel.Omega(:,:,modes)/wavemodel.f0));

RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);
RedundantCoefficients2D = shift(RedundantCoefficients);

if normalizeByGMAmplitude == 1
    wavemodel.InitializeWithGMSpectrum(1.0,'maxK', max(k)*2/3);
    A = shift( abs(wavemodel.Amp_plus(:,:,modes)) );
end

% Create a structure to log the most resonant triad for each wavenumber
nothing = NaN(size(K));
resonance = struct('ii2',nothing,'jj2',nothing,'kk2',nothing,'k2',nothing,'l2',nothing,'m2',nothing,'omega2',nothing,'ii3',nothing,'jj3',nothing,'kk3',nothing,'k3',nothing,'l3',nothing,'m3',nothing,'omega3',nothing,'Amp',nothing);

% individual wave I want to find triads for...
iK = 70;
iL = 17;

for iK=1:length(k)
    for iL=1:length(l)        
        for iM=1:length(modes)
%             if RedundantCoefficients2D(iK,iL)==1
%                 continue;
%             end
            
            k_i = K(iK,iL,iM);
            l_i = L(iK,iL,iM);
            m_i = M(iK,iL,iM);
            omega_i = Omega(iK,iL,iM);
            
            % These are all possible combinations (including self and redundants) of
            % resultant vectors.
            k_p = k_i + K;
            l_p = l_i + L;
            m_p = m_i + M;
            omega_p = omega_i + Omega;
            
            % negative omega??!?
            omega_m = abs(omega_i - Omega);
            
            % let's make these column vectors
            k_p = reshape(k_p,[],1);
            l_p = reshape(l_p,[],1);
            m_p = reshape(m_p,[],1);
            omega_p = reshape(omega_p,[],1);
            omega_m = reshape(omega_m,[],1);
            
            % first, find the location of the resultant vector (might be nan!)
            k_bin = discretize(k_p+dk/2,k);
            l_bin = discretize(l_p+dl/2,l);
            m_bin = discretize(m_p+dm/2,m);
            
            % find the resultant vectors that are *not* nan...
            omega_bin = nan(size(k_bin));
            validIndices = find(~isnan(k_bin) & ~isnan(l_bin) & ~isnan(m_bin));
            linearIndices = sub2ind(size(Omega),k_bin(validIndices),l_bin(validIndices),m_bin(validIndices));
            %...and find the associated frequencies
            omega_bin(validIndices) = Omega(linearIndices);
            
            if normalizeByGMAmplitude == 1
                A2 = reshape(A,[],1);
                A3 = nan(size(k_bin));
                A3(validIndices) = A(linearIndices);
                % maybe not multiplied by A3?
                Amp = max(A2.*A3./abs(omega_bin-omega_p),A2.*A3./abs(omega_bin-omega_m));
            else
                Amp = max(1./abs(omega_bin-omega_p), 1./abs(omega_bin-omega_m) );
            end
            Amp(logical(reshape(RedundantCoefficients2D,[],1))) = NaN;
            [maxAmp,maxAmpIndex] = max(Amp,[],'omitnan');
            if isnan(maxAmp)
                continue;
            end
            
            [ii,jj,kk] = ind2sub(size(K),maxAmpIndex);
            resonance.ii2(iK,iL,iM) = ii;
            resonance.jj2(iK,iL,iM) = jj;
            resonance.kk2(iK,iL,iM) = kk;
            resonance.k2(iK,iL,iM) = K(ii,jj,kk);
            resonance.l2(iK,iL,iM) = L(ii,jj,kk);
            resonance.m2(iK,iL,iM) = M(ii,jj,kk);
            resonance.omega2(iK,iL,iM) = Omega(ii,jj,kk);
            
            resonance.ii3(iK,iL,iM) = k_bin(maxAmpIndex);
            resonance.jj3(iK,iL,iM) = l_bin(maxAmpIndex);
            resonance.kk3(iK,iL,iM) = m_bin(maxAmpIndex);
            resonance.k3(iK,iL,iM) = K(k_bin(maxAmpIndex),l_bin(maxAmpIndex),m_bin(maxAmpIndex));
            resonance.l3(iK,iL,iM) = L(k_bin(maxAmpIndex),l_bin(maxAmpIndex),m_bin(maxAmpIndex));
            resonance.m3(iK,iL,iM) = M(k_bin(maxAmpIndex),l_bin(maxAmpIndex),m_bin(maxAmpIndex));
            resonance.omega3(iK,iL,iM) = omega_bin(maxAmpIndex);
            
            resonance.omega12(iK,iL,iM) = omega_p(maxAmpIndex);
            
            resonance.Amp(iK,iL,iM) = maxAmp;
        end
    end
end

[sortedResonance,sortOrder] = sort(reshape(resonance.Amp,[],1),'descend');
originalIndex = (1:length(sortedResonance))';
originalIndex = originalIndex(sortOrder);

firstNonNan = find(~isnan(sortedResonance),1,'first');

for iResonant=(firstNonNan+(0:10))
    maxAmp = sortedResonance(iResonant);
    maxAmpIndex = originalIndex(iResonant);
    
    [ii,jj,kk] = ind2sub(size(K),maxAmpIndex);
    %     fprintf('\nResonant pair %d is (i,j,k)=(%d,%d,%d) added to (i,j,k)=(%d,%d,%d) to get (i,j,k)=(%d,%d,%d)\n',iResonant-firstNonNan+1,ii,jj,kk,resonance.ii2(ii,jj,kk),resonance.jj2(ii,jj,kk),resonance.kk2(ii,jj,kk),resonance.ii3(ii,jj,kk),resonance.jj3(ii,jj,kk),resonance.kk3(ii,jj,kk));
    i0=length(k)/2+1;
    j0=length(l)/2+1;
    fprintf('\nResonant pair %d is (i,j,k)=(%d,%d,%d) added to (i,j,k)=(%d,%d,%d) to get (i,j,k)=(%d,%d,%d)\n',iResonant-firstNonNan+1,ii-i0,jj-j0,kk,resonance.ii2(ii,jj,kk)-i0,resonance.jj2(ii,jj,kk)-j0,resonance.kk2(ii,jj,kk),resonance.ii3(ii,jj,kk)-i0,resonance.jj3(ii,jj,kk)-j0,resonance.kk3(ii,jj,kk));
    fprintf('This is (k,l,omega): (%.2g, %.2g,%.2f f0) + (%.2g, %.2g,%.2f f0) = (%.2g, %.2g,%.2f f0)\n',K(ii,jj,kk),L(ii,jj,kk),Omega(ii,jj,kk),resonance.k2(ii,jj,kk),resonance.l2(ii,jj,kk),resonance.omega2(ii,jj,kk),resonance.k3(ii,jj,kk),resonance.l3(ii,jj,kk),resonance.omega3(ii,jj,kk));
end

figure
jpcolor(k,l,log10(resonance.Amp(:,:,1)'));
% hold on
% scatter(k_i+dk/2,l_i+dl/2,25,0*[1 1 1],'filled')

return



figure
jpcolor(k,l,RedundantCoefficients2D(:,:,1).')
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

