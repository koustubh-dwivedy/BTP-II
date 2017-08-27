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