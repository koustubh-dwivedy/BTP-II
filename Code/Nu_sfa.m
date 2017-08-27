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