function [m_used, r, v, u, m] = lander(t_f, r_0, v_0, r_N, v_N, p)
% constants and vehicle parameters
    N = 100;
	dt = t_f/(N-1);            		% period of calculation
	a = 1/(p.Isp*p.g0);			% alpha used in mass calculations
	gs = 4;						% glides slope constraint
	m_t = p.m_d + p.m_f;

	cvx_begin QUIET
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
				v(:,k+1) == v(:,k) + ((dt/2)*(u(:,k) + u(:,k+1))) +(p.g_plan*dt);
				z(1,k+1) == z(1,k) - (((a*dt)/2)*(s(1,k) + s(1,k+1)));
			end

			for k=1:N
				norm(u(:,k)) <= s(1,k);
				z_0 = m_t - (a*p.rho2*dt*(k-1));
				m_1 = p.rho1/z_0;
				m_2 = p.rho2/z_0;
				z1 = log(m_t-(a*p.rho1*dt*(k-1)));	
				z0 = log(z_0);
				z(1,k) >= z0;
				z(1,k) <= z1;
				s(1,k) <= m_2*(1 - (z(1,k) - z0));
				s(1,k) >= m_1*(1 - (z(1,k) - z0) + (((z(1,k) - z0)^2)/2));
			end
	cvx_end


	if strcmp(cvx_status, 'Solved')
	    m = exp(z);
	    m_used = m(1) - m(N);
	elseif strcmp(cvx_status, 'Infeasible')
	    m_used = m_t;
	else
	    fprintf('Error! %s', cvx_status);
    end

end






