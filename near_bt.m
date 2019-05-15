clear all; close all;

% dynamics near the Bogdanov-Takens bifurcation in the persistent sodium plus potassium model
% note: we assume tau_n(v)=1, C=1;

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

p=[-31.64,-31.55,-31.55,-31.70,-31.70;-79.42,-79.40,-79.44,-79.40,-79.44];
x_traj=[-56,-57,-58,-59.5];
figure
for i=1:4
    
%     vh_n = p(1,i);
%     E_l = p(2,i);
    % Compute nullclines and vector field
    %vVec=(-90:0.001:20); % vector spanning voltage values (v-dimension)
    %nVec=(0:0.001:0.7)'; % vector spanning n-gate values (n-dimension)
    vVec=(-70:0.01:-40); % vector spanning voltage values (v-dimension)
    nVec=(0.01:0.0001:0.3)'; % vector spanning n-gate values (n-dimension)
    [X,Y]=meshgrid(vVec,nVec); % create 2D grid of downsampled values
    v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
    n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
    v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
    n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline
    
    
%     sorted = vVec(sort(abs(v_null-n_null)))
%     fixed_inds = find(abs(v_null-n_null)<0.000003); % indices of fixed pts
    win = [-58.4,-58.3;-57.9,-57.6;-43,-42];
    for fixed_i=1:3
        win_inds=[find(vVec>win(fixed_i,1),1,'first'),find(vVec<win(fixed_i,2),1,'last')];
        [B,dist_inds] = sort(abs(v_null(win_inds)-n_null(win_inds)));
        fixed_inds(fixed_i) = win_inds(dist_inds(1));
    end
    
    
    % set trajectory starting points
    startx=x_traj(i)*ones(1,35);
    starty=(0.0272:0.0002:0.034)-0.004*(i-1); % iteratively decrease the starting y-coordinates so that they stay centered around the nullclines
    
    % plot nullclines and trajectories
    subplot(4,2,i*2-1)
    plot(vVec,v_null,'k:',vVec,n_null,'k-.')
    axis([-65 -50 nVec(1) 0.06])
    hold on
    scatter(vVec(fixed_inds),n_null(fixed_inds),'ro'); % equilibria
    streamline(X,Y,v_dot,n_dot,startx,starty) % trajectories
    ylimits=get(gca,'ylim');
    line([x_traj(i) x_traj(i)],ylimits,'Color','green')
    hold off
    legend('V-nullcline','n-nullcline')
    xlabel('membrane potential, V (mV)')
    ylabel('K^+ activation, n')
    
    subplot(4,2,i*2)
    plot(vVec,v_null,'k:',vVec,n_null,'k-.')
    hold on
    scatter(vVec(fixed_inds),n_null(fixed_inds),'ro'); % equilibria
    streamline(X,Y,v_dot,n_dot,startx,starty)
    hold off
    axis([vVec(find(vVec>=-61,1,'first')) vVec(find(vVec<=-55,1,'last')) 0.016 0.026])
    xlabel('membrane potential, V (mV)')
    ylabel('K^+ activation, n')
end

% % Simulation
% T=50; % total time
% dt=0.1; % simulation time-step (ms)
% tVec=(dt:dt:T); % time vector
% numTime=length(tVec); % number of timepoints
% 
% % inital values in phase-space
% v=-59.5*ones(1,numTime);
% n=0.02*ones(1,numTime);
% 
% % injected current
% I=I_base*ones(1,numTime); % baseline injected current
% p_w=round(5/dt); % pulse width
% I(round(10/dt):round(10/dt)+p_w)=-10;% pulse of input DC current
% 
% % forward Euler method
% for i=1:numTime-1
%     m_inf = 1/(1+exp((vh_m-v(i))/k_m));
%     n_inf = 1/(1+exp((vh_n-v(i))/k_n));
%     n(i+1) = n(i)+dt*(n_inf-n(i));
%     v(i+1) = v(i)+dt*(I(i)-g_l*(v(i)-E_l)-g_na*m_inf*(v(i)-E_na)-g_k*n(i)*(v(i)-E_k));
% end
% 
% % plot spiking results
% figure
% yyaxis left
% plot(tVec,v)
% ylabel('potential, V (mV)')
% axis([0 max(tVec) min(v)-0.15/0.80*(max(v)-min(v)) max(v)+0.05/0.80*(max(v)-min(v))])
% yyaxis right
% plot(tVec,I)
% axis([0 max(tVec) min(I)-0.05/0.05*(max(I)-min(I)) max(I)+0.90/0.05*(max(I)-min(I))])
% yticks(unique(I))
% ylabel('injected current, I ({\mu}A/cm^2)')
% xlabel('time (ms)')
