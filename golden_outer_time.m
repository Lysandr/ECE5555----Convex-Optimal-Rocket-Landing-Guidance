


function [tv, m_used, r, v, u, m] = golden_outer_time()
    tic;
    p.g0 = 9.807;				% Earth gravity, m/s^2
	p.g_plan = [0; 0;-3.711];	% Other planet gravity, m/s^2
	p.tv_max = 25;				% maximum TVC angle
	p.Isp = 255;				% specific impulse, s
	p.m_d = 1500;             	% drymass, kg
	p.m_f = 500;              	% fuel mass, kg
	p.Ft =  22000;        		% thrust
	p.rho2 = p.Ft;              	% thrust, Newtons
	p.rho1 = 0.15 * p.Ft;       	% lowest throttleability, Newtons
	r_0 = [-2000; 1500; 2000];	% position vector, m
	v_0 = [50; 70; -75];		% velocity vector, m/s
	r_N =[0; 0; 0];				% terminal position, m
	v_N =[0; 0; 0];				% terminal velocity, m
	t_fmin = 0;
	t_fmax = 200;

    obj_fun = @(t_f)(lander(t_f, r_0, v_0, r_N, v_N, p));
    options = optimset('TolX',0.5,'Display','iter');
    tf_opt = fminbnd(obj_fun, t_fmin, t_fmax, options);
    
    % Re-run optimal case
    [m_used, r, v, u, m] = GFOLD_fix_time(tf_opt, N, r0, v0, rf, vf, m_wet, theta, p);
    tv = linspace(0, tf_opt, N);
    plot_run3D(tv, r, v, u, m)
    toc;
end
