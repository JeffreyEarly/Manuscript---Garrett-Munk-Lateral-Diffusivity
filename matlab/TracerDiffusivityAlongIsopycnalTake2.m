% Notes on January 19th, 2020
%
% Following along an isopycnal does improve the measure of concentration.
% In fact, you find that the concentration of dye along an isopycnal varies
% in 1 part per 10000, rather than 1 part per 100 when you compute
% concentration along depth. So, that's a good fix! However, there is
% *nothing* else that gets rid of the oscillations in the second moment. It
% might be fair to say that these oscillations are simply real.

runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end

outputfile = sprintf('%s_tracer_column_iso.mat',file);
outputfile_ex = sprintf('%s_tracer_column_ex.mat',file);

if ~exist(outputfile,'file')
    WM = WintersModel(file);
    wavemodel = WM.wavemodel;
    
    dx = WM.wavemodel.x(2)-WM.wavemodel.x(1);
    dy = WM.wavemodel.y(2)-WM.wavemodel.y(1);
    dz = WM.wavemodel.z(2)-WM.wavemodel.z(1);
    X = WM.wavemodel.X;
    Y = WM.wavemodel.Y;
    Z = WM.wavemodel.Z;

%     dV = dx*dy*dz*ones(1,1,length(wavemodel.z));
%     dV(1) = dV(1)/2;
%     dV(end) = dV(end)/2;
%     


dX = zeros(size(X));
dY = zeros(size(X));
dZ = zeros(size(X));

% background density
RHOBAR = zeros(size(X)) + shiftdim(wavemodel.rhobar,-2);
    
    nFiles = WM.NumberOf3DOutputFiles;
    fileIncrements = 1:nFiles;
    
    nT = length(fileIncrements);
    t = zeros(nT,1);
    m = zeros(nT,wavemodel.Nz);
    m_x = zeros(nT,wavemodel.Nz);
    m_y = zeros(nT,wavemodel.Nz);
    m_xx = zeros(nT,wavemodel.Nz);
    m_xy = zeros(nT,wavemodel.Nz);
    m_yy = zeros(nT,wavemodel.Nz);
    
    iIso = 22;
    % dye concentration along an isopycnal
    Cn = zeros(nT,wavemodel.Nx, wavemodel.Ny);
    XC = X(:,:,iIso);
    YC = Y(:,:,iIso);
    XisoC = zeros(nT,wavemodel.Nx, wavemodel.Ny);
    YisoC = zeros(nT,wavemodel.Nx, wavemodel.Ny);
    
    startTime = datetime('now');
    for iFile = 1:length(fileIncrements)
        if mod(iFile,10) == 2
            timePerStep = (datetime('now')-startTime)/(iFile-1);
            timeRemaining = (nT-iFile+1)*timePerStep;
            fprintf('reading values time step %d of %d to file. Estimated finish time %s (%s from now)\n', iFile, nT, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
        end
        
        file = WM.PathOf3DOutputFileAtIndex(fileIncrements(iFile));
        [rho_prime,s2,t(iFile)] = WM.VariableFieldsFrom3DOutputFileAtIndex(iFile,'rho_prime','s2','time');
        
        % the density on an X,Y,Z grid
        rho = shiftdim(wavemodel.rhobar,-2) + rho_prime;
        
        % density forced to be monotonic. Hopefully won't have changed too much.
        RHO_monotonic = sort(rho,3,'descend');
        
        Z_isopycnal = zeros(size(X));
        for iX=1:size(X,1)
            for iY=1:size(X,2)
                x = reshape(RHO_monotonic(iX,iY,:),[],1);
                v = reshape(Z(iX,iY,:),[],1);
                xq = reshape(RHOBAR(iX,iY,:),[],1);
                Z_isopycnal(iX,iY,:) = interp1(x,v,xq,'linear');
            end
        end
        
        % distance between isopycnals
        dZ(:,:,1) = Z_isopycnal(:,:,2)-Z_isopycnal(:,:,1);
        dZ(:,:,2:(end-1)) = (Z_isopycnal(:,:,3:end)-Z_isopycnal(:,:,1:(end-2)))/2;
        dZ(:,:,end) = Z_isopycnal(:,:,end)-Z_isopycnal(:,:,end-1);
        
        % distance along isopycnals
        Xiso = cumsum(sqrt(dx^2 + dZ.^2),1)-dx;
        Yiso = cumsum(sqrt(dy^2 + dZ.^2),2)-dy;
        
        % distance increment, ignoring end points
        dX(2:(end-1),:,:) = (Xiso(3:end,:,:) - Xiso(1:(end-2),:,:))/2;
        dY(:,2:(end-1),:) = (Yiso(:,3:end,:) - Yiso(:,1:(end-2),:))/2;
        
        dV = dX.*dY.*dZ;
        
        s2_isopycnal = interpn(X,Y,Z,s2,X,Y,Z_isopycnal,'linear');
        
        % std(m,0,1)./mean(m,1) = O(1e-4), so pretty damn good. Without
        % handling isopycnals, we get O(1e-2)
        m(iFile,:) = squeeze(sum(sum(s2_isopycnal.*dV,1),2));
        
        
        Cn(iFile,:,:) = s2_isopycnal(:,:,iIso).*dV(:,:,iIso)/m(iFile,iIso);
        XisoC(iFile,:,:) = Xiso(:,:,iIso);
        YisoC(iFile,:,:) = Yiso(:,:,iIso);
        
        % Need to compute the length along the isopycnal?
        % sum(abs(Z_isopycnal(:,1,64)-wavemodel.z(64))) is 2 km, which is
        % the size of the oscillations.
        m_x(iFile,:) = shiftdim(sum(sum(X.*s2_isopycnal.*dV,1),2),1)./m(iFile,:);
        m_y(iFile,:) = shiftdim(sum(sum(Y.*s2_isopycnal.*dV,1),2),1)./m(iFile,:);
        
        Xc = X-shiftdim(m_x(iFile,:),-1);
        Yc = Y-shiftdim(m_y(iFile,:),-1);

        m_xx(iFile,:) = shiftdim(sum(sum(Xc.*Xc.*s2_isopycnal.*dV,1),2),1)./m(iFile,:);
        m_xy(iFile,:) = shiftdim(sum(sum(Xc.*Yc.*s2_isopycnal.*dV,1),2),1)./m(iFile,:);
        m_yy(iFile,:) = shiftdim(sum(sum(Yc.*Yc.*s2_isopycnal.*dV,1),2),1)./m(iFile,:);

    end

    save(outputfile_ex,'t','Cn','XC', 'YC','XisoC','YisoC');
    save(outputfile,'t','m','m_x', 'm_y','m_xx','m_xy','m_yy');
else
    load(outputfile);
end


D2 = (m_xx+m_yy)/4;

% D2 = m_yy/2;
kappa = zeros(size(D2,2),1);
kappa_err = zeros(size(D2,2),1);
for iZ=1:size(D2,2)
    y = D2(:,iZ);
    [D2_coeff,D2_err] = LinearLeastSquaresFit(t,y);
    kappa(iZ) = D2_coeff(2);
    kappa_err(iZ) = D2_err(2);
end

x = [kappa+3*kappa_err; flip(kappa-3*kappa_err)];
y = [wavemodel.z; flip(wavemodel.z)];

figure('Name','LateralDiffusivityOfTracerVsDepth')
fill(x,y,[1 1 1]*0.8), hold on
plot(kappa,wavemodel.z,'LineWidth',2);
xlabel('lateral diffusivity (m^2/s)')
ylabel('depth (m)')
ylim([min(wavemodel.z) max(wavemodel.z)])
xlim([0 max(x)*1.05])
