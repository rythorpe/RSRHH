clear all; close all;

% persistent sodium plus potassium model
% note: we assume tau_n(v)=1, C=1;

%% Integrator: postinhibitory spike
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

T=50; % total time
dt=0.1; % simulation time-step (ms)
tVec=(dt:dt:T); % time vector
numTime=length(tVec); % number of timepoints

% inital values in phase-space
v=-59.5*ones(1,numTime);
n=0.02*ones(1,numTime);

% injected current
I=4.3*ones(1,numTime); % baseline injected current
p_w=round(5/dt); % pulse width
I(round(10/dt):round(10/dt)+p_w)=-10;% pulse of input DC current

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
axis([0 max(tVec) min(v)-3/10*(max(v)-min(v)) max(v)+1/10*(max(v)-min(v))])
yyaxis right
plot(tVec,I)
axis([0 max(tVec) min(I)-1/2*(max(I)-min(I)) max(I)+8*(max(I)-min(I))])
yticks(unique(I))
ylabel('injected current, I ({\mu}A/cm^2)')
xlabel('time (ms)')


%% Resonator
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
vVec=(-90:1:20); % vector spanning voltage values (v-dimension)
nVec=(0:0.01:0.7)'; % vector spanning n-gate values (n-dimension)
[X,Y]=meshgrid(downsample(vVec,2),downsample(nVec,2)); % create 2D grid of downsampled values
v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline

% Plot nullclines and vector field
figure
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
axis([vVec(1) vVec(end) 0 0.7])
hold on
h=quiver(X,Y,v_dot,n_dot,'AutoScaleFactor',2); % downsample vVec and plot vector field 
hs = get(h,'MaxHeadSize');
set(h,'MaxHeadSize',hs/1000)
startx=[-80,-60,-60,-60,-35];
starty=[0.04,0.01,0.06,0.11,0.04];
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
legend('V-nullcline','n-nullcline')
xlabel('membrane potential, V (mV)')
ylabel('K^+ activation, n')
axes('Position',[.58 .2 .18 .18])
box on
plot(vVec,v_null,'k:',vVec,n_null,'k-.')
hold on
streamline(X,Y,v_dot,n_dot,startx,starty)
hold off
axis([vVec(find(vVec>=-65,1,'first')) vVec(find(vVec<=-55,1,'last')) 0.01 0.05])
box off

% Simulation
T=50; % total time
dt=0.1; % simulation time-step (ms)
tVec=(dt:dt:T); % time vector
numTime=length(tVec); % number of timepoints

% inital values in phase-space
v=-59.5*ones(1,numTime);
n=0.02*ones(1,numTime);

% injected current
I=I_base*ones(1,numTime); % baseline injected current
p_w=round(5/dt); % pulse width
I(round(10/dt):round(10/dt)+p_w)=-10;% pulse of input DC current

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
axis([0 max(tVec) min(v)-3/10*(max(v)-min(v)) max(v)+1/10*(max(v)-min(v))])
yyaxis right
plot(tVec,I)
axis([0 max(tVec) min(I)-1/2*(max(I)-min(I)) max(I)+8*(max(I)-min(I))])
yticks(unique(I))
ylabel('injected current, I ({\mu}A/cm^2)')
xlabel('time (ms)')
