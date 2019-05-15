clear all; close all;

% persistent sodium plus potassium model
% note: we assume tau_n(v)=1, C=1;

% simulation parameters
T=300; % total time (ms)
dt=0.1; % simulation time-step (ms)
tVec=(dt:dt:T); % time vector
numTime=length(tVec); % number of timepoints
pw=1; % pulse width (ms)
ipi=[3.8 7.6 11.4]; % inter-pulse-interval (ms)
pulse_inds=[(find(tVec>=10,1,'first'):find(tVec<=10+pw,1,'last')),...
    (find(tVec>=10+ipi(1),1,'first'):find(tVec<=10+ipi(1)+pw,1,'last')),...
    (find(tVec>=110,1,'first'):find(tVec<=110+pw,1,'last')),...
    (find(tVec>=110+ipi(2),1,'first'):find(tVec<=110+ipi(2)+pw,1,'last')),...
    (find(tVec>=210,1,'first'):find(tVec<=210+pw,1,'last')),...
    (find(tVec>=210+ipi(3),1,'first'):find(tVec<=210+ipi(3)+pw,1,'last'))];

%%% integrator %%%
% parameters
g_l=8; % conductance constant (leak)
g_na=20;
g_k=10;
E_l=-79.42; % Nernst potential (leak)
E_na=60;
E_k=-90;
vh_m=-20; % steady-state half-life potential (m-gate)
vh_n=-31;
k_m=15; % sigma function slope factor (m-gate)
k_n=7;
I_base=4.3; % baseline injected current

% Compute nullclines and vector field
vVec=(-90:0.1:20); % vector spanning voltage values (v-dimension)
nVec=(0:0.001:0.7)'; % vector spanning n-gate values (n-dimension)
[X,Y]=meshgrid(downsample(vVec,50),downsample(nVec,50)); % create 2D grid of downsampled values
v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline

% Plot nullclines and vector field
figure
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
axis([vVec(1) vVec(end) 0 0.7])
hold on
scatter([-40.5,-57.5,-59.4],[0.2047,0.02219,0.017],'ro'); % equilibria
h=quiver(X,Y,v_dot,n_dot,'AutoScaleFactor',2); % downsample vVec and plot vector field 
hs = get(h,'MaxHeadSize');
set(h,'MaxHeadSize',hs/1000)
startx=[-60,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5];
starty=[0.04,0.6,0.6,0.6,0.6,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.01,0.01,0.01,0.01,];
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
title('Integrator neuron')
legend('V-nullcline','n-nullcline')
xlabel('membrane potential, V (mV)')
ylabel('K^+ activation, n')
axes('Position',[.58 .2 .18 .18])
box on
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
hold on
scatter([-40.5,-57.5,-59.4],[0.2047,0.02219,0.017],'ro'); % equilibria
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
axis([vVec(find(vVec>=-61,1,'first')) vVec(find(vVec<=-54,1,'last')) 0.016 0.024])
box off

% inital values in phase-space
v=-59.4*ones(1,numTime);
n=0.017*ones(1,numTime);

% injected current
I=I_base*ones(1,numTime); % baseline injected current
I(pulse_inds)=I_base+0.68;% % add pulse current

% forward Euler method
for i=1:numTime-1
    m_inf = 1/(1+exp((vh_m-v(i))/k_m));
    n_inf = 1/(1+exp((vh_n-v(i))/k_n));
    n(i+1) = n(i)+dt*(n_inf-n(i));
    v(i+1)=v(i)+dt*(I(i)-g_l*(v(i)-E_l)-g_na*m_inf*(v(i)-E_na)-g_k*n(i)*(v(i)-E_k));
end

% plot spiking results
figure
yyaxis left
plot(tVec,v)
ylabel('potential, V (mV)')
axis([0 max(tVec) min(v)-0.15/0.80*(max(v)-min(v)) max(v)+0.05/0.80*(max(v)-min(v))])
yyaxis right
plot(tVec,I)
axis([0 max(tVec) min(I)-0.05/0.05*(max(I)-min(I)) max(I)+0.90/0.05*(max(I)-min(I))])
yticks(unique(I))
title('Integrator neuron')
ylabel('injected current, I ({\mu}A/cm^2)')
xlabel('time (ms)')


%%% resonator %%%
% parameters
g_l=8; % conductance constant (leak)
g_na=20;
g_k=10;
E_l=-79.42; % Nernst potential (leak)
E_na=60;
E_k=-90;
vh_m=-20; % steady-state half-life potential (m-gate)
vh_n=-34;
k_m=15; % sigma function slope factor (m-gate)
k_n=7;
I_base=7; % baseline injected current

% Compute nullclines and vector field
vVec=(-90:0.1:20); % vector spanning voltage values (v-dimension)
nVec=(0:0.001:0.7)'; % vector spanning n-gate values (n-dimension)
[X,Y]=meshgrid(downsample(vVec,50),downsample(nVec,50)); % create 2D grid of downsampled values
v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline

% Plot nullclines and vector field
figure
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
axis([vVec(1) vVec(end) 0 0.7])
hold on
scatter(-59.4,0.02587,'ro'); % equilibria
h=quiver(X,Y,v_dot,n_dot,'AutoScaleFactor',2); % downsample vVec and plot vector field 
hs = get(h,'MaxHeadSize');
set(h,'MaxHeadSize',hs/1000)
startx=[-60,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5];
starty=[0.04,0.6,0.6,0.6,0.6,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.01,0.01,0.01,0.01,];
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
title('Resonator neuron')
legend('V-nullcline','n-nullcline')
xlabel('membrane potential, V (mV)')
ylabel('K^+ activation, n')
axes('Position',[.58 .2 .18 .18])
box on
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
hold on
scatter(-59.4,0.02587,'ro'); % equilibria
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
axis([vVec(find(vVec>=-65,1,'first')) vVec(find(vVec<=-55,1,'last')) 0.02 0.04])
box off

% inital values in phase-space
v=-59.4*ones(1,numTime);
n=0.02587*ones(1,numTime);

% injected current
I=I_base*ones(1,numTime); % baseline injected current
I(pulse_inds)=I_base+1.9;% add pulse current

% forward Euler method
for i=1:numTime-1
    m_inf = 1/(1+exp((vh_m-v(i))/k_m));
    n_inf = 1/(1+exp((vh_n-v(i))/k_n));
    n(i+1) = n(i)+dt*(n_inf-n(i));
    v(i+1) = v(i)+dt*(I(i)-g_l*(v(i)-E_l)-g_na*m_inf*(v(i)-E_na)-g_k*n(i)*(v(i)-E_k));
end

% plot spiking results
figure
yyaxis left
plot(tVec,v)
ylabel('potential, V (mV)')
axis([0 max(tVec) min(v)-0.15/0.80*(max(v)-min(v)) max(v)+0.05/0.80*(max(v)-min(v))])
yyaxis right
plot(tVec,I)
axis([0 max(tVec) min(I)-0.05/0.05*(max(I)-min(I)) max(I)+0.90/0.05*(max(I)-min(I))])
yticks(unique(I))
title('Resonator neuron')
ylabel('injected current, I ({\mu}A/cm^2)')
xlabel('time (ms)')


%%% BT %%%
% parameters
g_l=8; % conductance constant (leak)
g_na=20;
g_k=10;
E_l=-79.42; % Nernst potential (leak)
E_na=60;
E_k=-90;
vh_m=-20; % steady-state half-life potential (m-gate)
vh_n=-31.64;
k_m=15; % sigma function slope factor (m-gate)
k_n=7;
I_base=5; % baseline injected current

% Compute nullclines and vector field
vVec=(-90:0.1:20); % vector spanning voltage values (v-dimension)
nVec=(0:0.001:0.7)'; % vector spanning n-gate values (n-dimension)
[X,Y]=meshgrid(downsample(vVec,50),downsample(nVec,50)); % create 2D grid of downsampled values
v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline

% Plot nullclines and vector field
figure
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
axis([vVec(1) vVec(end) 0 0.7])
hold on
scatter([-42.3,-57.7,-58.4],[0.179,0.02359,0.0214],'ro'); % equilibria
h=quiver(X,Y,v_dot,n_dot,'AutoScaleFactor',2); % downsample vVec and plot vector field 
hs = get(h,'MaxHeadSize');
set(h,'MaxHeadSize',hs/1000)
startx=[-60,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5,-85,-55,-25,5];
starty=[0.04,0.6,0.6,0.6,0.6,0.4,0.4,0.4,0.4,...
    0.2,0.2,0.2,0.2,0.01,0.01,0.01,0.01,];
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
title('Transition (near BT) neuron')
legend('V-nullcline','n-nullcline')
xlabel('membrane potential, V (mV)')
ylabel('K^+ activation, n')
axes('Position',[.58 .2 .18 .18])
box on
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
hold on
scatter([-42.3,-57.7,-58.4],[0.179,0.02359,0.0214],'ro'); % equilibria
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
axis([vVec(find(vVec>=-65,1,'first')) vVec(find(vVec<=-55,1,'last')) 0.02 0.04])
box off

% inital values in phase-space
v=-58.396*ones(1,numTime);
n=0.0214*ones(1,numTime);

% injected current
I=I_base*ones(1,numTime); % baseline injected current
I((find(tVec>=10,1,'first'):find(tVec<=10+pw,1,'last')))=I_base+0.8;% add pulse current

% forward Euler method
for i=1:numTime-1
    m_inf = 1/(1+exp((vh_m-v(i))/k_m));
    n_inf = 1/(1+exp((vh_n-v(i))/k_n));
    n(i+1) = n(i)+dt*(n_inf-n(i));
    v(i+1) = v(i)+dt*(I(i)-g_l*(v(i)-E_l)-g_na*m_inf*(v(i)-E_na)-g_k*n(i)*(v(i)-E_k));
end

% plot spiking results
figure
yyaxis left
plot(tVec,v)
ylabel('potential, V (mV)')
axis([0 max(tVec) min(v)-0.15/0.80*(max(v)-min(v)) max(v)+0.05/0.80*(max(v)-min(v))])
yyaxis right
plot(tVec,I)
axis([0 max(tVec) min(I)-0.05/0.05*(max(I)-min(I)) max(I)+0.90/0.05*(max(I)-min(I))])
yticks(unique(I))
title('Transition (near BT) neuron')
ylabel('injected current, I ({\mu}A/cm^2)')
xlabel('time (ms)')
