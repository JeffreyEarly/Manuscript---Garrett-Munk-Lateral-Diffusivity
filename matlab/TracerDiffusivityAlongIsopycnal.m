runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_11/';
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

x = wavemodel.x;
y = wavemodel.y;
z = wavemodel.z;
dx = WM.wavemodel.x(2)-WM.wavemodel.x(1);
dy = WM.wavemodel.y(2)-WM.wavemodel.y(1);
dz = WM.wavemodel.z(2)-WM.wavemodel.z(1);
X = WM.wavemodel.X;
Y = WM.wavemodel.Y;
Z = WM.wavemodel.Z;

dV = dx*dy*dz;

[t,u,v,w,rho_prime,s2] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime','s2');
m = sum(sum(sum(s2)))*dV;
m_x = sum(sum(sum(X.*s2)))*dV/m;
m_y = sum(sum(sum(Y.*s2)))*dV/m;
m_z = sum(sum(sum(Z.*s2)))*dV/m;

m_i = find(wavemodel.x<m_x,1,'last');
m_j = find(wavemodel.y<m_y,1,'last');
m_k = find(wavemodel.z<m_z,1,'last');

rho = shiftdim(wavemodel.rhobar,-2) + rho_prime;

isopycnal = rho(m_i,m_j,m_k);

rho_xz = squeeze(rho(:,m_j,:));
tracer_xz = squeeze(s2(:,m_j,:));

C = contourc(wavemodel.x,wavemodel.z,rho_xz.',isopycnal*[1 1]);
numVertices = C(2,1);
x = C(1,1+(1:numVertices)).';
z = C(2,1+(1:numVertices)).';

figure
subplot(1,2,1)
pcolor(wavemodel.x,wavemodel.z,rho_xz.'),shading flat
hold on,
plot(x,z,'k','LineWidth',2)
subplot(1,2,2)
pcolor(wavemodel.x,wavemodel.z,tracer_xz.'),shading flat
hold on,
plot(x,z,'k','LineWidth',2)

s = [0; cumsum(sqrt(diff(x).^2 + diff(z).^2))];