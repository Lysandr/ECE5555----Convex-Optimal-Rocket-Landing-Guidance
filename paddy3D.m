%% Padraig Lysandrou
% Implementation of the Becet, Ploen Algorithm
clear all;
close all;

% constants and vehicle parameters
g0 = 9.807;					% Earth gravity, m/s^2
g_plan = [0; 0;-3.711];		% Other planet gravity, m/s^2
tv_max = 25;				% maximum TVC angle
Isp = 255;					% specific impulse, s
m_d = 1500;             	% drymass, kg
m_f = 500;              	% fuel mass, kg
m_t = m_f + m_d;			% total mass, kg
Ft =  22000;        		% thrust
rho2 = Ft;              	% thrust, Newtons
rho1 = 0.15 * Ft;       	% lowest throttleability, Newtons
% initial and final conditions -- dynamics  [x y z]
r_0 = [-2000; 1500; 2000];	% position vector, m
v_0 = [50; 70; -75];		% velocity vector, m/s
r_N =[0; 0; 0];				% terminal position, m
v_N =[0; 0; 0];				% terminal velocity, m
% other timing and constraint garbage
t_f = 100;					% final time horizon (not necessarily optimal)
dt = 1;             		% period of calculation
N = 1+(t_f/dt);				% calculation steps
a = 1/(Isp*g0);				% alpha used in mass calculations
gs = 4;						% glides slope constraint


% the discretized problem 4
% maximize the terminal mass
cvx_solver SEDUMI
cvx_begin
	variables u(3,N) z(1,N) s(1,N) r(3,N) v(3,N)
	minimize(-z(N))			% objective function
	subject to 				% constraints and dynamics
		r(:,1) == r_0;		% position IC
		v(:,1) == v_0;		% velocity IC
		r(:,N) == r_N; 		% position TC
		v(:,N) == v_N;		% velocity TC
		z(1) == log(m_t);	% mass IC

		r(3,:) >= 0;							% plz don't crash into the ground
		%u(2,:) >= s.*cos(degtorad(tv_max));		% thrust vector control constraint
		%r(1,:) <= r(2,:)/tan(degtorad(gs));		% glide slope
        z(:) >= 0;

		for  k = 1:N-1
			r(:,k+1) == r(:,k) + ((dt/2)*(v(:,k) + v(:,k+1))) +(((dt^2)/12)*(u(:,k+1) - u(:,k)));
			v(:,k+1) == v(:,k) + ((dt/2)*(u(:,k) + u(:,k+1))) +(g_plan*dt);
			z(1,k+1) == z(1,k) - (((a*dt)/2)*(s(1,k) + s(1,k+1)));
		end

		for k=1:N
			norm(u(:,k)) <= s(1,k);
			z_0 = m_t - (a*rho2*dt*(k-1));
			m_1 = rho1/z_0;
			m_2 = rho2/z_0;
			z1 = log(m_t-(a*rho1*dt*(k-1)));	
			z0 = log(z_0);
			z(1,k) >= z0;
			z(1,k) <= z1;
			s(1,k) <= m_2*(1 - (z(1,k) - z0));
			s(1,k) >= m_1*(1 - (z(1,k) - z0) + (((z(1,k) - z0)^2)/2));
		end

cvx_end

%% plotski
plot_all_data3D(dt,t_f,r,v,u,z)







