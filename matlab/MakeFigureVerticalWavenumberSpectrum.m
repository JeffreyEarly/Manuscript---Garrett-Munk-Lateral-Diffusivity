scaleFactor = 1;
LoadFigureDefaults;

runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2019_01/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end

WM = WintersModel(file);
wavemodel = WM.wavemodel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%

[t,u,v,zeta] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','zeta');

N0 = WM.wavemodel.internalModes.N0;

b = 1300;

% Let's be extra formal here...
xi = cumtrapz(wavemodel.z,sqrt(wavemodel.N2));
Lxi = xi(end)-xi(1);
s = (b/Lxi)*xi;
Ls = s(end)-s(1);
% although with constant stratification 1300m deep, this should be the same as z.

% now let's imagine we only sampled the top half of the domain
zIndices = floor(wavemodel.Nz/2):wavemodel.Nz;
% zIndices = 1:wavemodel.Nz; % try the whole water column to start with.
z = wavemodel.z(zIndices);
Lz = max(z)-min(z);
s_p = s(zIndices);
Ls_p = max(s_p)-min(s_p);

u_p = u(:,:,zIndices);
v_p = v(:,:,zIndices);

% [u_bar] = m/s
[u_bar,~] = CosineTransformForward(s_p,u_p,3);
[v_bar,m] = CosineTransformForward(s_p,v_p,3);
dm = m(2)-m(1);
j = 2*b*m;

% [E] = m^3/s^2
E = (1.3e3)*(1.3e3)*(1.3e3)*(5.2e-3)*(5.2e-3)*(6.3e-5);

% [S_u] m^3/s^2
S_u = Ls_p*squeeze(mean(mean(u_bar.*conj(u_bar) + v_bar.*conj(v_bar),2),1));

% variance sanity check
Eu_dz = mean(mean(trapz(z,u_p.*u_p + v_p.*v_p,3),2),1)/Lz;
Eu_dm = (S_u(1)/2 + sum(S_u(2:end)))*dm;

% full depth should be...
Eu_gm = (3/2)*E/b;

S_u_wkb = S_u;


Hj = S_u_wkb*df/Eu_gm;

j_star = 3;
H1 = (j_star+(1:3000)).^(-5/2);
H_norm = 1/sum(H1);
H = @(j) H_norm*(j_star + j).^(-5/2);

j_theory = 1:max(j);
m_theory = j_theory/(2*b);
figure
scatter(m,Hj), hold on
scatter(m_theory,H(j_theory))
xlog, ylog

return
figure
scatter(j(2:end),Hj(2:end)./H(j(2:end)))