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
function f = T_BB(height, percent_cloud_cover)
	if height >= 6000
		f = 204.4;
	else
		f = 214.4 - 0.20*percent_cloud_cover;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vector = F(t,y)
	g = 9.81;
	m_sys = 3.0054;
	m_helium = 0.409;
	m_pf = 0.815;
	m_sf = 0.157;
	m_freon = 1.000;
	liq_r114_density = 1518.093; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane
	R_helium = 2077;
	R_freon = 48.644;
	c_D = 0.8;
	c_M = 0.5;
	L = 3; % Length of Secondary Balloon %Approximated from image of ALICE0/D in paper
	W = 0.9*3/3.1; % Width of Secondary Balloon %Approximated from image of ALICE0/D in paper
	T_boil = 276.6;
	c_p_he = 5193; %J/(Kg.K)
	c_p_polyethylene = 1650; %http://www.tainstruments.com/pdf/literature/TA227.pdf
	%c_p_polyethylene = ; % Assumed Mol weight of % Came by averaging Cp of branched polyethylene (LDPE) at 256K(260) and 223K(220) (US Standard Atmosphere Temp at 5km and 10km resp). There's scope of improving this as a constant value is incorrect.
	c_p_r114_gas = 689.56; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane %Gas property at 15 degree Celsius
	latent_heat_vap_r114 = 135939; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane %Gas property at 15 degree Celsius
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
	liq_r114_vol = y(8)/liq_r114_density;
	vol_primary = temp_amb*m_helium*R_helium/pressure_amb;
	vol_secondary = temp_amb*(m_freon - y(8))*R_freon/pressure_amb;
	radius_primary = nthroot((3*vol_primary)/(4*pi),1/3);
	radius_secondary = nthroot((3*vol_secondary)/(4*pi),1/3);
	S_p = 4*pi*radius_primary*radius_primary;
	S_s = 4*pi*radius_secondary*radius_secondary;
	S_sg = 2*L*W;
	S_sl = (pi/2)*W*nthroot((8*liq_r114_vol)/(pi*W),1/2);

	percent_cloud_cover = ;

	G = 1396;
	sigma = 5.670367*(10^(-8));
	r_e = 0.18 + 0.0039*percent_cloud_cover;
	T_BB_val = T_BB(y(2), percent_cloud_cover);
	epsilon_pint = ;
	epsilon_pweff = ;
	epsilon_pgeff = ;
	epsilon_sint = ;
	epsilon_sgeff = ;
	epsilon_sleff = ;
	epsilon_sweff = ;
	CH_pgf = ;
	CH_pfa = ;
	CH_lfa = ;
	CH_sgf = ;
	alpha_pgeff = ;
	alpha_pweff = ;
	alpha_sgeff = ;
	alpha_sleff = ;
	alpha_sweff = ;

	q_dot_pg = (G*alpha_pgeff*(1 + r_e) + epsilon_pint*sigma*(y(3)^4 - y(4)^4) - CH_pgf*(y(3) - y(4)) - epsilon_pgeff*sigma*(y(3)^4) + epsilon_pgeff*sigma*(T_BB_val^4))*S_p;
	q_dot_pf = (G*alpha_pweff*(0.25 + 0.5*r_e) + epsilon_pint*sigma*(y(3)^4 - y(4)^4) + CH_pgf*(y(3) - y(4)) + CH_pfa*(temp_amb - y(4)) - epsilon_pweff*sigma*(y(4)^4) + epsilon_pweff*sigma*(T_BB_val^4))*S_p;
	q_dot_sg = (G*alpha_sgeff*(1 + r_e) + epsilon_sint*sigma*(y(5)^4 - y(7)^4) - epsilon_sgeff*sigma*(y(5)^4) + epsilon_sgeff*sigma*(T_BB_val^4))*S_s - CH_sgf*(y(5) - y(7))*S_sg;
	q_dot_sl = (G*alpha_sleff*(0.5 + r_e) + epsilon_sleff*sigma*(T_BB_val^4 - y(6)^4) + CH_lfa*(temp_amb - y(6)))*S_sl;
	q_dot_sf = (G*alpha_sweff*(0.25 + 0.5*r_e) + epsilon_sint*sigma*(y(5)^4 - y(7)^4) - epsilon_sweff*sigma*(y(7)^4) + epsilon_sweff*sigma*(T_BB_val^4))*S_s + (CH_sgf*(y(5) - y(7)) + CH_sfa*(temp_amb - y(7)))*S_sg;

	vector = [
	(g*(density_amb*(vol_primary + vol_secondary) - m_sys) - 0.5*density_amb*c_D*y(1)*abs(y(1))*pi*radius_primary*radius_primary)/(m_sys + c_M*density_amb*(vol_primary + vol_secondary));
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