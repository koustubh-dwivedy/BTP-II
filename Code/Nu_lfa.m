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