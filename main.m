global state;
state = [0, 0, 1]; %because T_boil =276.6K and initial conditions are at 300K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = switch_master(T_sl, T_sg, T_boil, m_sl, m_s)
	global state;
	m_buffer = 0.05*m_s;
	if (state(1) == 0 && state(2) == 1 && state(3) == 0)
		if (T_sl > T_boil)
			state = [1, 0, 1];
		end
	elseif (state(1) == 1 && state(2) == 0 && state(3) == 1)
		if (m_sl < m_buffer)
			state = [0, 0, 1];
		elseif (m_sl > (m_s - m_buffer))
			state = [0, 1, 0];
		end
	elseif (state(1) == 0 && state(2) == 0 && state(3) == 1)
		if (T_sg < T_boil)
			state = [1, 0, 1];
		end
	end
	f = [state(1), state(2), state(3)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = density(x)
	[a b c d] = atmosisa(x);
	f = d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = pressure(x)
	[a b c d] = atmosisa(x);
	f = c;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = temperature(x)
	[a b c d] = atmosisa(x);
	f = a;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vector = F(t,y)
	g = 9.81;
	m_sys = 3.0054;
	m_helium = 0.409;
	R_helium = 2077;
	m_freon = 1.000;
	R_freon = 48.644;
	c_D = 0.8;
	c_M = 0.5;
	T_boil = 276.6;
	% y(1)-> dz/dt
	% y(2)-> z
	% y(3)-> T_pg
	% y(4)-> T_pf
	% y(5)-> T_sg
	% y(6)-> T_sl
	% y(7)-> T_sf
	% y(8)-> m_sl
	[a, b, c] = switch_master(y(6), y(5), T_boil, y(8), m_freon);
	vector = [
	(g*(density(y(2))*temperature(y(2))*(m_helium*R_helium + (m_freon - y(8))*R_freon)/pressure(y(2)) - m_sys) - 0.5*density(y(2))*c_D*y(1)*abs(y(1))*pi*nthroot((3*m_helium*R_helium*temperature(y(2)))/(4*pi*pressure(y(2))),2/3))/(m_sys + c_M*density(y(2))*temperature(y(2))*(m_helium*R_helium + (m_freon - y(8))*R_freon)/pressure(y(2)));
	;
	];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%