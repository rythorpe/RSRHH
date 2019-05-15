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

x_traj=[-56,-57,-58,-59.5]; % x-coordinates for initial trajectories in each plot
figure
for i=1:4
    % Compute nullclines and vector field
    vVec=(-70:0.01:-40); % vector spanning voltage values (v-dimension)
    nVec=(0.01:0.0001:0.3)'; % vector spanning n-gate values (n-dimension)
    [X,Y]=meshgrid(vVec,nVec); % create 2D grid of downsampled values
    v_dot = I_base-g_l.*(X-E_l)-g_na./(1+exp((vh_m-X)/k_m)).*(X-E_na)-g_k*Y.*(X-E_k); % dv/dt
    n_dot = 1./(1+exp((vh_n-X)/k_n))-Y; % dn/dt
    v_null = (I_base-g_l*(vVec-E_l)-g_na./(1+exp((vh_m-vVec)/k_m)).*(vVec-E_na))./(g_k*(vVec-E_k)); % v-nullcline
    n_null = 1./(1+exp((vh_n-vVec)/k_n)); % n-nullcline
    
    % find fixed points
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
    
    % plot nullclines and trajectories zoomed in near lower equilibria
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
