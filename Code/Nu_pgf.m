function f = Nu_pgf(rho_he, eta_he, g, dia, T_film, T_gas, Pr_he)
	Gr = ((rho_he/eta_he)^2)*g*(dia^3)*(abs(T_film - T_gas)/T_gas);
	Rn = Gr*Pr_he;
	if Rn < (1.34681*(10^8))
		f = 2.5*(2 + 0.6*(Rn^0.25));
	else
		f = 0.325*(Rn^(1/3));
	end
end