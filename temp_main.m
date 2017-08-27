global state;
state = [0, 0, 1]; %because T_boil =276.6K and initial conditions are at 300K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a b c] = switch_master(T_sl, T_sg, T_boil, m_sl, m_s)
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
	a = state(1);
	b = state(2);
	c = state(3);
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
function f = cloud_cover(height)
	f = 68;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Nu_pfa(V_p, rho_air, dia, vel, eta_air, g, T_film, T_gas, Pr_a)
	Re = rho_air*dia*vel/eta_air;
	Gr = ((rho_air/eta_air)^2)*g*(dia^3)*(abs(T_film - T_gas)/T_gas);
	if V_p < 53800
		Nu_1 = 0.37*(Re^0.6);
	else
		Nu_1 = 0.74*(Re^0.6);
	end
	Nu_2 = 2 + 0.6*((Gr*Pr_a)^0.25);
	f = max([Nu_1, Nu_2]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Nu_pgf(rho_he, eta_he, g, dia, T_film, T_gas, Pr_he)
	Gr = ((rho_he/eta_he)^2)*g*(dia^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_he;
	if Rn < (1.34681*(10^8))
		f = 2.5*(2 + 0.6*(Rn^0.25));
	else
		f = 0.325*(Rn^(1/3));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Nu_lfa(rho_air, dia, vel, eta_air, g, T_film, T_gas, Pr_a)
	Re = rho_air*dia*vel/eta_air;
	Gr = ((rho_air/eta_air)^2)*g*(dia^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_a;
	if Re < 1297.742
		Nu_1 = 0.38 + 0.44*(Re^0.5);
	else
		Nu_1 = 0.22*(Re^0.6);
	end
	Nu_2 = (0.6 + 0.322*(Rn^(1/6)))^2;
	f = max([Nu_1, Nu_2]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Nu_sgf(rho_r114, eta_r114, g, L, T_film, T_gas, Pr_r114)
	Gr = ((rho_r114/eta_r114)^2)*g*(L^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_r114;
	Nu_T = 0.515*(Rn^0.25);
	Nu_l = 2.8/(log(1 + 2.8/Nu_T));
	Nu_t = 0.103*(Rn^(1/3));
	f = (Nu_t^6 + Nu_l^6)^(1/6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Nu_sfa(rho_air, eta_air, g, len, T_film, T_gas, vel, Pr_a)
	Re = rho_air*len*vel/eta_air;
	Gr = ((rho_air/eta_air)^2)*g*(len^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_a;
	if Re < 487508.3
		Nu_1 = 0.5924*(Re^0.5);
	else
		Nu_1 = 0.033*(Re^0.8) - 758.3;
	end

	Nu_T = 0.515*(Rn^0.25);
	Nu_l = 2.8/(log(1 + 2.8/Nu_T));
	Nu_t = 0.103*(Rn^(1/3));
	Nu_2 = (Nu_t^6 + Nu_l^6)^(1/6);
	f = max([Nu_1, Nu_2]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vector = F(t,y)
	liq_r114_density = 1518.093; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane
	g = 9.81;
	R_helium = 2077;
	R_freon = 48.644;
	c_D = 0.8;
	c_M = 0.5;
	L = 3; % Length of Secondary Balloon %Approximated from image of ALICE0/D in paper
	W = 0.9*3/3.1; % Width of Secondary Balloon %Approximated from image of ALICE0/D in paper
	T_boil = 276.6;
	latent_heat_vap_r114 = 135939; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane %Gas property at 15 degree Celsius

	eta_he = 0.0000187; %http://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html
	eta_air = 0.0000173; %http://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html
	eta_r114 = 0.000010412; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane

	Pr_a = 0.71;
	Pr_he = 0.68; %http://www.mhtl.uwaterloo.ca/old/onlinetools/airprop/airprop.html
	Pr_r114 = 0.739; %http://www.ethermo.us/ShowDetail41.htm

	m_sys = 3.0054;
	m_helium = 0.409;
	m_pf = 0.815;
	m_sf = 0.157;
	m_freon = 1.000;

	c_p_he = 5193; %J/(Kg.K)
	c_p_polyethylene = 1650; %http://www.tainstruments.com/pdf/literature/TA227.pdf
	%c_p_polyethylene = ; % Assumed Mol weight of % Came by averaging Cp of branched polyethylene (LDPE) at 256K(260) and 223K(220) (US Standard Atmosphere Temp at 5km and 10km resp). There's scope of improving this as a constant value is incorrect.
	c_p_r114_gas = 689.56; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane %Gas property at 15 degree Celsius
	c_p_r114_liq = 982; %Cp exp (experimental) at 20.5deg Celsius. http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=910720

	k_air = 0.022; %http://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
	k_he = 0.142; %http://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
	k_r114 = 0.009702; %https://encyclopedia.airliquide.com/12-dichloro-1122-tetrafluoroethane

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

	rho_he = m_helium/vol_primary;
	rho_r114 = (m_freon - y(8))/vol_secondary;

	radius_primary = nthroot((3*vol_primary)/(4*pi),1/3);
	radius_secondary = nthroot((3*vol_secondary)/(4*pi),1/3);

	S_p = 4*pi*radius_primary*radius_primary;
	S_s = 4*pi*radius_secondary*radius_secondary;
	S_sg = 2*L*W;
	S_sl = (pi/2)*W*nthroot((8*liq_r114_vol)/(pi*W),1/2);

	percent_cloud_cover = cloud_cover(y(2));

	G = 1396;
	sigma = 5.670367*(10^(-8));
	r_e = 0.18 + 0.0039*percent_cloud_cover;
	T_BB_val = T_BB(y(2), percent_cloud_cover);

	r_pw = 0.127;
	r_sw = 0.127;
	r_pwsol = 0.114;
	r_swsol = 0.114;

	tau_pw = 0.98; %at wavenumber 500 (10000/18.3) (ref actual paper and this->) http://people.csail.mit.edu/jaffer/FreeSnell/polyethylene.html
	tau_sw = 0.98; %at wavenumber 500 (10000/18.3) (ref actual paper and this->) http://people.csail.mit.edu/jaffer/FreeSnell/polyethylene.html
	tau_pwsol = 0.82; %at wavenumber 4000 (10000000/2500) (ref actual paper and this->) http://people.csail.mit.edu/jaffer/FreeSnell/polyethylene.html
	tau_swsol = 0.82; %at wavenumber 4000 (10000000/2500) (ref actual paper and this->) http://people.csail.mit.edu/jaffer/FreeSnell/polyethylene.html

	epsilon_pg = 0.0003;
	epsilon_pw = 0.031;
	epsilon_sg = 0.60;
	epsilon_sw = 0.031;
	epsilon_sl = 0.95;

	alpha_pg = 0.003;
	alpha_pw = 0.001;
	alpha_sg = 0.00184;
	alpha_sw = 0.001;
	alpha_sl = 0.003;

	epsilon_pint = epsilon_pg*epsilon_pw/(1 - r_pw*(1 - epsilon_pg));
	epsilon_pweff = epsilon_pw*(1 + (tau_pw*(1 - epsilon_pg))/(1 - r_pw*(1 - epsilon_pg)));
	epsilon_pgeff = epsilon_pg*tau_pw/(1 - r_pw*(1 - epsilon_pg));
	epsilon_sint = epsilon_sg*epsilon_sw/(1 - r_sw*(1 - epsilon_sg));
	epsilon_sweff = epsilon_sw*(1 + (tau_sw*(1 - epsilon_sg))/(1 - r_sw*(1 - epsilon_sg)));
	epsilon_sgeff = epsilon_sg*tau_sw/(1 - r_sw*(1 - epsilon_sg));
	epsilon_sleff = epsilon_sl*tau_sw/(1 - r_sw*(1 - epsilon_sl)); % MIGHT BE WRONG. SIMPLY EXTRAPOLATED

	CH_pgf = Nu_pgf(rho_he, eta_he, g, 2*radius_primary, y(4), y(3), Pr_he)*k_he/(2*radius_primary);
	CH_pfa = Nu_pfa(vol_primary, density_amb, 2*radius_primary, abs(y(1)), eta_air, g, y(4), temp_amb, Pr_a)*k_air/(2*radius_primary);
	CH_lfa = Nu_lfa(density_amb, (S_sl*2/(pi*W)), abs(y(1)), eta_air, g, y(7), temp_amb, Pr_a)*k_air/(S_sl*2/(pi*W));
	CH_sgf = Nu_sgf(rho_r114, eta_r114, g, L, y(7), y(5), Pr_r114)*k_r114/L;
	CH_sfa = Nu_sfa(density_amb, eta_air, g, L, y(7), temp_amb, abs(y(1)), Pr_a)*k_air/L;

	alpha_pgeff = alpha_pg*tau_pwsol/(1 - r_pwsol*(1 - alpha_pg));
	alpha_pweff = alpha_pw*(1 + (tau_pwsol*(1 - alpha_pg))/(1 - r_pwsol*(1 - alpha_pg)));
	alpha_sgeff = alpha_sg*tau_swsol/(1 - r_swsol*(1 - alpha_sg));
	alpha_sweff = alpha_sw*(1 + (tau_swsol*(1 - alpha_sg))/(1 - r_swsol*(1 - alpha_sg)));
	alpha_sleff = alpha_sl*tau_swsol/(1 - r_swsol*(1 - alpha_sl)); % MIGHT BE WRONG. SIMPLY EXTRAPOLATED

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
[t y]=ode45(F,[0 7200], [0 500 300 300 300 300 300 0]);
% y(1)-> dz/dt
% y(2)-> z
% y(3)-> T_pg
% y(4)-> T_pf
% y(5)-> T_sg
% y(6)-> T_sl
% y(7)-> T_sf
% y(8)-> m_sl