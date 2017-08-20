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
	m_pf = 0.815;
	m_sf = 0.157;
	R_helium = 2077;
	m_freon = 1.000;
	R_freon = 48.644;
	c_D = 0.8;
	c_M = 0.5;
	T_boil = 276.6;
	c_p_he = 5193; %J/(Kg.K)
	c_p_polyethylene = 1650; %http://www.tainstruments.com/pdf/literature/TA227.pdf
	%c_p_polyethylene = ; % Assumed Mol weight of % Came by averaging Cp of branched polyethylene (LDPE) at 256K(260) and 223K(220) (US Standard Atmosphere Temp at 5km and 10km resp). There's scope of improving this as a constant value is incorrect.
	c_p_r114_gas = 689.56; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane Gas property at 15 degree Celsius
	latent_heat_vap_r114 = 135939; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane Gas property at 15 degree Celsius
	c_p_r114_liq = 982; %Cp exp (experimental) at 20.5deg Celsius. http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=910720
	% y(1)-> dz/dt
	% y(2)-> z
	% y(3)-> T_pg
	% y(4)-> T_pf
	% y(5)-> T_sg % multiply by c
	% y(6)-> T_sl % multiply by b
	% y(7)-> T_sf
	% y(8)-> m_sl % multiply by a
	[a, b, c] = switch_master(y(6), y(5), T_boil, y(8), m_freon);
	density_amb = density(y(2));
	temp_amb = temperature(y(2));
	pressure_amb = pressure(y(2));
	vol_primary = temp_amb*m_helium*R_helium/pressure_amb;
	vol_secondary = temp_amb*(m_freon - y(8))*R_freon/pressure_amb;
	q_dot_pg = ;
	q_dot_pf = ;
	q_dot_sg = ;
	q_dot_sl = ;
	q_dot_sf = ;

	vector = [
	(g*(density_amb*(vol_primary + vol_secondary) - m_sys) - 0.5*density_amb*c_D*y(1)*abs(y(1))*pi*nthroot((3*vol_primary)/(4*pi),2/3))/(m_sys + c_M*density_amb*(vol_primary + vol_secondary));
	y(1);
	(q_dot_pg - density_amb*g*y(1)*vol_primary)/(m_helium*c_p_he);
	q_dot_pf/(c_p_polyethylene*m_pf);
	c*(q_dot_sg - density_amb*g*y(1)*vol_secondary)/(c_p_r114_gas*((m_freon - y(8))));
	b*q_dot_sl/(y(8)*c_p_r114_liq);
	q_dot_sf/(m_sf*c_p_polyethylene);
	a*q_dot_sl/latent_heat_vap_r114];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t y]=ode45(F,[0 7200], []);