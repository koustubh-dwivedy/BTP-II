function f = Nu_sgf(rho_r114, eta_r114, g, L, T_film, T_gas, Pr_r114)
	Gr = ((rho_r114/eta_r114)^2)*g*(L^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_r114;
	Nu_T = 0.515*(Rn^0.25);
	Nu_l = 2.8/(log(1 + 2.8/Nu_T));
	Nu_t = 0.103*(Rn^(1/3));
	f = (Nu_t^6 + Nu_l^6)^(1/6);
end