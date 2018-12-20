runtype = 'linear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
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

nFiles = WM.NumberOf3DOutputFiles;
[t0,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t0,u,v,rho_prime);

U = max(max(max(abs(sqrt(u.^2 + v.^2)))));
T = min(abs(2*pi./wavemodel.Omega(:)));


[t_p,x_p,y_p,z_p] = WM.ParticleTrajectories;

iParticle = 101;
p0 = [x_p(1,iParticle), y_p(1,iParticle), z_p(1,iParticle)-wavemodel.Lz];

f = @(t,y) wavemodel.DrifterVelocityAtTimePositionVector(t,y,'exact');

% Let's do fixed step size integrator.
cfl = 0.25;
advectiveDT = cfl*(wavemodel.x(2)-wavemodel.x(1))/U;
oscillatoryDT = T/8;
if advectiveDT < oscillatoryDT
    fprintf('Using the advective dt: %.2f\n',advectiveDT);
    deltaT = advectiveDT;
else
    fprintf('Using the oscillatory dt: %.2f\n',oscillatoryDT);
    deltaT = oscillatoryDT;
end

t_indices = 1:50;
t_in = t_p(t_indices);

% p2 = ode4(f,t_in, p0');
[t,p2] = ode45(f,t_in, p0,odeset('RelTol',1e-5,'AbsTol',1e-8));

x2 = p2(:,1);
y2 = p2(:,2);
z2 = p2(:,3);

figure
subplot(3,1,1)
plot(t_in,x2)
hold on, plot(t_in,x_p(t_indices,iParticle))
subplot(3,1,2)
plot(t_in,y2)
hold on, plot(t_in,y_p(t_indices,iParticle))
subplot(3,1,3)
plot(t_in,z2)
hold on, plot(t_in,z_p(t_indices,iParticle)-wavemodel.Lz)