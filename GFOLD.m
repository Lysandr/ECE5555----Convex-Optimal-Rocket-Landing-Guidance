% Neal Bhasin
% 2015-04-21
% G-FOLD outer time-optimization routine. 
%
% Inputs:
% N     : Number of knot points in discrete optimization problem
% r0    : Initial position [m]   |  rf : Final position [m]
% v0    : Initial velocity [m/s] |  vf : Final velocity [m/s]
% m_wet : Initial total mass [kg]
% theta : Thrust pointing limit from vertical [deg]
% p : Vehicle/planet parameters
% p.m_dry        : Vehicle mass without fuel [kg]
% p.Isp          : Specific impulse [s]
% p.g            : Planet gravity vector [m/s^2]
% p.max_throttle : Max open throttle [.%]
% p.min_throttle : Min open throttle [.%]
% p.T_max        : Max total thrust force at 1.0 throttle [N]
% p.phi          : The cant angle of thrusters [deg]
%
% Outputs:
% tv     : Vector (1xN) of time at knot points [s]
% m_used : Fuel mass used [kg]
% r      : Position(t) [m]
% v      : Velocity(t) [m/s]
% u      : Command Acceleration(t) [m/s^2]
% m      : Mass(t) [kg]

function [tv, m_used, r, v, u, m] = GFOLD(N, r0, v0, rf, vf, m_wet, theta, p)
    tic;
    g0 = 9.80665; % Standard earth gravity [m/s^2]
    alpha = 1 / (p.Isp * g0 * cosd(p.phi));
    r1 = p.min_throttle * p.T_max * cosd(p.phi);
    r2 = p.max_throttle * p.T_max * cosd(p.phi);
    tf_min = p.m_dry * norm(vf - v0) / r2;
    tf_max = (m_wet - p.m_dry) / (alpha * r1);
    
    cvx_solver SEDUMI
    obj_fun = @(t)( GFOLD_fix_time(t, N, r0, v0, rf, vf, m_wet, theta, p) );
    options = optimset('TolX',0.5,'Display','iter');
    tf_opt = fminbnd(obj_fun, tf_min, tf_max, options);
    
    % Re-run optimal case
    [m_used, r, v, u, m] = GFOLD_fix_time(tf_opt, N, r0, v0, rf, vf, m_wet, theta, p);
    tv = linspace(0, tf_opt, N);
    plot_run2D(tv, r, v, u, m)
    toc;
end
