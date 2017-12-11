% Neal Bhasin
% 2015-04-21
% G-FOLD inner fuel-optimization routine.

% See GFOLD.m for description of arguments.
function [m_used, r, v, u, m] = GFOLD_fix_time(tf, N, r0, v0, rf, vf, m_wet, theta, p)
dt = tf / (N - 1);

g0 = 9.80665; % Standard earth gravity [m/s^2]
alpha = 1 / (p.Isp * g0 * cosd(p.phi));
r1 = p.min_throttle * p.T_max * cosd(p.phi);
r2 = p.max_throttle * p.T_max * cosd(p.phi);

cvx_begin QUIET
    % Parameterize trajectory position, velocity, thrust acceleration, ln mass
    variables r(2,N) v(2,N) u(2,N) z(1,N) s(1,N)
    % Maximize ln of final mass -> Minimize fuel used
    maximize(z(N))
    subject to
        % Initial condition constraints
        r(:,1) == r0;
        v(:,1) == v0;
        z(1) == log(m_wet);
        % Terminal condition constraints
        r(:,N) == rf;
        v(:,N) == vf;
        % Dynamical constraints
        for i=1:N-1
            % Position / Velocity
            v(:,i+1) == v(:,i) + dt*p.g + (dt/2)*(u(:,i) + u(:,i+1));
            r(:,i+1) == r(:,i) + (dt/2)*(v(:,i) + v(:,i+1)) + ...
                (dt^2/12)*(u(:,i+1) - u(:,i));
            % Mass
            z(i+1) == z(i) - (alpha*dt/2)*(s(i) + s(i+1));
        end
        % Thrust limit, mass flow limit
        for i=1:N
            norm(u(:,i)) <= s(i);
            % Feasible/conservative Taylor series expansion
            z0_term = m_wet - alpha * r2 * (i-1) * dt;
            z1_term = m_wet - alpha * r1 * (i-1) * dt;
            z0 = log(z0_term);
            z1 = log(z1_term);
            mu_1 = r1 / z0_term;
            mu_2 = r2 / z0_term;
            % Quadratic lower bound, linear upper bound
            s(i) >= mu_1 * (1 - (z(i) - z0) + (1/2)*(z(i) - z0)^2);
            s(i) <= mu_2 * (1 - (z(i) - z0));
            % Impose physical extremal mass limits
            z(i) >= z0;
            z(i) <= z1;
        end
        % Thrust pointing constraint
        u(2,:) >= s .* cosd(theta);
        % No sub-surface flight
        r(2,1:N-1) >= 0;
cvx_end
if strcmp(cvx_status, 'Solved')
    m = exp(z);
    m_used = m(1) - m(N);
elseif strcmp(cvx_status, 'Infeasible')
    m_used = m_wet;
else
    fprintf('Error! %s', cvx_status);
end
end
